###
# EPPS 7390 Bayesian Final Project
# Brennan Stout
# 12/10/2025
###

# 1. SETUP & DATA PREPARATION --------------------------------------------------

library(tidyverse)
library(lubridate)
library(rstanarm)
library(mvtnorm)
library(gridExtra)
library(dplyr)
if (requireNamespace("INLA", quietly = TRUE)) library(INLA)

# Load Data
data <- read.csv("https://raw.githubusercontent.com/BrennanStout/epps7390/refs/heads/main/WHO_Measles_12_25.csv")

# Preprocess
df_model <- data %>%
  group_by(year, month) %>%
  summarise(measles_total = sum(measles_total, na.rm = TRUE), .groups = 'drop') %>%
  arrange(year, month) %>%
  mutate(
    date = make_date(year, month, 1),
    time = 1:n(),
    log_lag1 = lag(log(measles_total + 1)),
    sin12 = sin(2 * pi * time / 12),
    cos12 = cos(2 * pi * time / 12),
    sin6 = sin(4 * pi * time / 12),
    cos6 = cos(4 * pi * time / 12),
    time_idx = time,
    month_idx = month
  ) %>%
  filter(!is.na(log_lag1))

# Container for all predictions
results_df <- data.frame()

########################### 
# MODEL 1: MCMC (rstanarm)

fit_mcmc <- stan_glm(
  measles_total ~ time + log_lag1 + sin12 + cos12 + sin6 + cos6,
  data = df_model,
  family = neg_binomial_2(link = "log"),
  chains = 2, iter = 2000, seed = 42, refresh = 0
)


# Extract Predictions
pred_mcmc <- posterior_predict(fit_mcmc)
system.time(
mcmc_res <- df_model %>%
  dplyr::select(date, measles_total) %>%
  mutate(
    Model = "MCMC",
    Mean_Est = colMeans(pred_mcmc),
    Lower_CI = apply(pred_mcmc, 2, quantile, probs = 0.025),
    Upper_CI = apply(pred_mcmc, 2, quantile, probs = 0.975)
  ))
results_df <- rbind(results_df, mcmc_res)

print(summary(fit_mcmc))
prior_summary(fit_mcmc)

# 5. Visualization
# Generate posterior predictions
# posterior_linpred returns the linear predictor (log scale), we want the expected count (response scale)
# posterior_epred returns the expectation (mean) of the response distribution
post_means <- posterior_epred(fit_mcmc)

# Calculate credible intervals for the mean
df_model$fit_mean <- colMeans(post_means)
df_model$fit_lower <- apply(post_means, 2, quantile, probs = 0.025)
df_model$fit_upper <- apply(post_means, 2, quantile, probs = 0.975)

# Plot
ggplot(df_model, aes(x = date)) +
  geom_point(aes(y = measles_total), alpha = 0.5, size = 2) +
  geom_line(aes(y = fit_mean), color = "firebrick", size = 1) +
  geom_ribbon(aes(ymin = fit_lower, ymax = fit_upper), fill = "firebrick", alpha = 0.2) +
  labs(
    title = "Measles Estimation: MCMC",
    y = "Total Cases",
    x = "Date"
  ) +
  theme_minimal()
# ==============================================================================
# MODEL 2: INLA
# ==============================================================================
cat("\n--- Running Model 2: INLA ---\n")

system.time(
if (requireNamespace("INLA", quietly = TRUE)) {
  formula_inla <- measles_total ~ 1 + f(time_idx, model = "rw2") + f(month_idx, model = "seasonal", season.length = 12)
  res_inla <- inla(formula_inla, family = "nbinomial", data = df_model,
                   control.predictor = list(compute = TRUE, link = 1))
  
  inla_res <- df_model %>%
    dplyr::select(date, measles_total) %>%
    mutate(
      Model = "INLA",
      Mean_Est = res_inla$summary.fitted.values$mean,
      Lower_CI = res_inla$summary.fitted.values$`0.025quant`,
      Upper_CI = res_inla$summary.fitted.values$`0.975quant`
    )
  results_df <- rbind(results_df, inla_res)
} else {
  cat("Skipping INLA (package not found).\n")
})

ggplot(inla_res, aes(x = date)) +
  geom_point(aes(y = measles_total), color = "black", alpha = 0.5, size = 2) +
  geom_line(aes(y = Mean_Est), color = "purple", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +
  labs(
    title = "Measles Estimation: INLA",
    y = "Total Cases",
    x = "Date",
    subtitle = "Model: Negative Binomial ~ RW2(time) + Seasonal(month)"
  ) +
  theme_minimal()

# ==============================================================================
# MODEL 3: ABC (Rejection Sampling)
# ==============================================================================
cat("\n--- Running Model 3: ABC ---\n")

y_obs <- df_model$measles_total
n_obs <- length(y_obs)
n_draws_abc <- 50000; keep_n_abc <- 200

# Priors
priors_abc <- list(
  b0 = runif(n_draws_abc, 0, 3), b_time = runif(n_draws_abc, -0.01, 0.01),
  b_ar = runif(n_draws_abc, 0.7, 0.99), alpha = runif(n_draws_abc, 0.01, 0.5),
  b_s12 = runif(n_draws_abc, -0.5, 0.5), b_c12 = runif(n_draws_abc, -0.5, 0.5),
  b_s6 = runif(n_draws_abc, -0.3, 0.3), b_c6 = runif(n_draws_abc, -0.3, 0.3)
)

# Simulation
sim_mat <- matrix(0, n_draws_abc, n_obs)
sim_mat[,1] <- y_obs[1]
curr <- rep(y_obs[1], n_draws_abc)
# Vectors
v_s12 <- df_model$sin12; v_c12 <- df_model$cos12; v_s6 <- df_model$sin6; v_c6 <- df_model$cos6; v_t <- df_model$time

set.seed(123)
system.time(
for(t in 2:n_obs){
  log_mu <- priors_abc$b0 + priors_abc$b_time * v_t[t] + 
    priors_abc$b_s12 * v_s12[t] + priors_abc$b_c12 * v_c12[t] +
    priors_abc$b_s6 * v_s6[t] + priors_abc$b_c6 * v_c6[t] + 
    priors_abc$b_ar * log(curr + 1)
  mu <- exp(pmin(pmax(log_mu, -5), 15))
  curr <- rnbinom(n_draws_abc, size = 1/priors_abc$alpha, mu = mu)
  sim_mat[,t] <- curr
})

# Rejection
dists <- rowSums(abs(sweep(sim_mat, 2, y_obs, "-")))
best_idx <- order(dists)[1:keep_n_abc]
best_sims <- sim_mat[best_idx, ]

abc_res <- df_model %>%
  dplyr::select(date, measles_total) %>%
  mutate(
    Model = "ABC",
    Mean_Est = colMeans(best_sims),
    Lower_CI = apply(best_sims, 2, quantile, probs = 0.025),
    Upper_CI = apply(best_sims, 2, quantile, probs = 0.975)
  )
results_df <- rbind(results_df, abc_res)

ggplot(abc_res, aes(x = date)) +
  geom_point(aes(y = measles_total), color = "black", alpha = 0.5, size = 1.5) +
  geom_line(aes(y = Mean_Est), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "blue", alpha = 0.2) +
  labs(
    title = "Measles Estimation: ABC",
    subtitle = "Negative Binomial + AR(1) + Dual Seasonality",
    y = "Total Cases",
    x = "Date"
  ) +
  theme_minimal()
# ==============================================================================
# MODEL 4: BSL (Bayesian Synthetic Likelihood)
# ==============================================================================
cat("\n--- Running Model 4: BSL ---\n")

# Stats & Model
get_stats <- function(y) {
  c(mean(y), log(var(y)+1e-6), ifelse(length(y)>1, acf(y, plot=F)$acf[2], 0), 
    cor(y, v_t), cor(y, v_s12), cor(y, v_c12), cor(y, v_s6), cor(y, v_c6))
}
S_obs <- get_stats(y_obs)

sim_bsl <- function(p, n, init) {
  y <- numeric(n); y[1] <- init; cur <- init
  fix <- p[1] + p[2]*v_t + p[3]*v_s12 + p[4]*v_c12 + p[5]*v_s6 + p[6]*v_c6
  sz <- 1/(p[8]+1e-9)
  for(t in 2:n){
    mu <- exp(pmin(pmax(fix[t] + p[7]*log(cur+1), -5), 12))
    cur <- rnbinom(1, size=sz, mu=mu); y[t] <- cur
  }
  y
}

# MCMC
n_bsl <- 2000; cur_p <- c(1, 0, 0, 0, 0, 0, 0.8, 0.1)
cov_p <- diag(c(0.05, 0.001, 0.05, 0.05, 0.05, 0.05, 0.02, 0.01)^2)
samps_bsl <- matrix(0, n_bsl, 8)

comp_ll <- function(p) {
  s <- t(replicate(30, get_stats(sim_bsl(p, n_obs, y_obs[1]))))
  dmvnorm(S_obs, colMeans(s), cov(s)+diag(1e-6,8), log=TRUE)
}
cur_ll <- comp_ll(cur_p)

system.time(
for(i in 1:n_bsl){
  prop <- rmvnorm(1, cur_p, cov_p)[1,]
  if(prop[8]>0 && prop[7]>=0 && prop[7]<=1.05){
    pll <- comp_ll(prop)
    if(is.finite(pll) && log(runif(1)) < (pll - cur_ll)) { cur_p <- prop; cur_ll <- pll }
  }
  samps_bsl[i,] <- cur_p
})

# Posterior Pred
burn <- 500; post_samps <- samps_bsl[(burn+1):n_bsl, ]
bsl_preds <- t(apply(post_samps[sample(nrow(post_samps), 200), ], 1, function(p) sim_bsl(p, n_obs, y_obs[1])))

bsl_res <- df_model %>%
  dplyr::select(date, measles_total) %>%
  mutate(
    Model = "BSL",
    Mean_Est = colMeans(bsl_preds),
    Lower_CI = apply(bsl_preds, 2, quantile, probs = 0.025),
    Upper_CI = apply(bsl_preds, 2, quantile, probs = 0.975)
  )
results_df <- rbind(results_df, bsl_res)

ggplot(bsl_res, aes(x = date)) +
  geom_point(aes(y = measles_total), alpha = 0.5) +
  geom_line(aes(y = Mean_Est), color = "darkred", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "red", alpha = 0.2) +
  labs(title = "Measles Estimation: BSL", y = "Cases") +
  theme_minimal()
# ==============================================================================
# MODEL 5: Variational Inference (ADVI)
# ==============================================================================
cat("\n--- Running Model 5: VI ---\n")

system.time(
fit_vi <- stan_glm(
  measles_total ~ time + log_lag1 + sin12 + cos12 + sin6 + cos6,
  data = df_model, family = neg_binomial_2(link = "log"),
  algorithm = "meanfield", iter = 10000, seed = 42, refresh = 0
))

pred_vi <- posterior_predict(fit_vi)
vi_res <- df_model %>%
  dplyr::select(date, measles_total) %>%
  mutate(
    Model = "VI",
    Mean_Est = colMeans(pred_vi),
    Lower_CI = apply(pred_vi, 2, quantile, probs = 0.025),
    Upper_CI = apply(pred_vi, 2, quantile, probs = 0.975)
  )
results_df <- rbind(results_df, vi_res)

ggplot(vi_res, aes(x = date)) +
  geom_point(aes(y = measles_total), color = "black", alpha = 0.5) +
  geom_line(aes(y = Mean_Est), color = "purple", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +
  labs(
    title = "Measles Estimation: Variational Inference (ADVI)",
    subtitle = "Approximation: Mean-Field Gaussian",
    y = "Total Cases",
    x = "Date"
  ) +
  theme_minimal()
# ==============================================================================
# VISUALIZATION & COMPARISON
# ==============================================================================
cat("\n--- Generating Visualizations ---\n")

# 1. Individual Model Plots
# plot_model <- function(df_sub) {
#   model_name <- unique(df_sub$Model)
#   ggplot(df_sub, aes(x = date)) +
#     geom_point(aes(y = measles_total), color = "black", alpha = 0.4, size = 1) +
#     geom_line(aes(y = Mean_Est), color = "blue", size = 0.8) +
#     geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "blue", alpha = 0.2) +
#     labs(title = paste("Model Fit:", model_name), y = "Cases", x = "Date") +
#     theme_minimal()
# }
# 
# plots_list <- lapply(split(results_df, results_df$Model), plot_model)


# 2. Compare Fit (Overlay)
p_compare <- ggplot(results_df, aes(x = date)) +
  geom_point(aes(y = measles_total), color = "black", alpha = 0.2, size = 1) +
  geom_line(aes(y = Mean_Est, color = Model), size = 0.8) +
  labs(title = "Model Comparison: Fit Overlay", y = "Cases", x = "Date") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_compare)

# 3. Performance Metrics (RMSE & MAE)
perf_metrics <- results_df %>%
  group_by(Model) %>%
  summarise(
    RMSE = sqrt(mean((measles_total - Mean_Est)^2)),
    MAE = mean(abs(measles_total - Mean_Est)),
    Coverage_95 = mean(measles_total >= Lower_CI & measles_total <= Upper_CI)
  ) %>%
  pivot_longer(cols = c(RMSE, MAE, Coverage_95), names_to = "Metric", values_to = "Value")

p_metrics <- ggplot(perf_metrics, aes(x = Model, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Model Performance Metrics", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_metrics)

print(perf_metrics %>% pivot_wider(names_from = Metric, values_from = Value))

#####################

plot_data <- data %>%
  mutate(date = make_date(year, month)) %>%
  group_by(region, date) %>%
  summarise(total_cases = sum(measles_total, na.rm = TRUE)) %>%
  ungroup()

# Create the Stacked Bar Chart
ggplot(plot_data, aes(x = date, y = total_cases, fill = region)) +
  geom_col(position = "stack", width = 30) + 
  labs(
    title = "Total Measles Cases by Region Over Time",
    x = "Date",
    y = "Total Measles Cases",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")