#Loo difference comparison#####
library(loo)


log_lik_array <- stan_out$draws("log_lik", format = "matrix")

loo_result  <- loo(log_lik_array)
waic_result <- waic(log_lik_array)


# saveRDS(loo_result,  file = "Mhet_LOO.rds")
#saveRDS(waic_result, file = "Mhet_WAIC.rds")

#saveRDS(loo_result,  file = "Mhom_LOO.rds")
#saveRDS(waic_result, file = "Mhom_WAIC.rds")


# Read saved model

loo1  <- readRDS("./Output/Mhet_LOO.rds")
waic1 <- readRDS("./Output/Mhet_WAIC.rds")

loo2  <- readRDS("./Output/Mhom_LOO.rds")
waic2 <- readRDS("./Output/Mhom_WAIC.rds")


# LOO comparison
loo_comp <- loo_compare(loo1, loo2)
print(loo_comp)

# Difference ModelA - ModelB
delta_elpd <-
  loo1$estimates["elpd_loo","Estimate"] -
  loo2$estimates["elpd_loo","Estimate"]

# Correct paired SE from loo_compare()
se_delta_elpd <- abs(loo_comp[2,"se_diff"])

cat("LOO difference =", delta_elpd, "\n")
cat("SE of LOO difference =", se_delta_elpd, "\n")


# WAIC comparison

delta_waic <-
  waic1$estimates["waic","Estimate"] -
  waic2$estimates["waic","Estimate"]

# Pointwise WAIC differences
waic_point_diff <-
  waic1$pointwise[, "elpd_waic"] -
  waic2$pointwise[, "elpd_waic"]

# paired SE
se_delta_waic <-
  2 * sqrt(length(waic_point_diff) * var(waic_point_diff))

cat("WAIC difference =", delta_waic, "\n")
cat("SE of WAIC difference =", se_delta_waic, "\n")



# comparison table

comparison_table <- data.frame(
  
  Metric = c("WAIC", "elpd_loo"),
  
  ModelA_Estimate = c(
    waic1$estimates["waic","Estimate"],
    loo1$estimates["elpd_loo","Estimate"]
  ),
  
  ModelA_SE = c(
    waic1$estimates["waic","SE"],
    loo1$estimates["elpd_loo","SE"]
  ),
  
  ModelB_Estimate = c(
    waic2$estimates["waic","Estimate"],
    loo2$estimates["elpd_loo","Estimate"]
  ),
  
  ModelB_SE = c(
    waic2$estimates["waic","SE"],
    loo2$estimates["elpd_loo","SE"]
  ),
  
  Difference = c(
    delta_waic,
    delta_elpd
  ),
  
  SE_Difference = c(
    se_delta_waic,
    se_delta_elpd
  )
  
)

print(comparison_table)



#MCMC diagnostics for two states#####

library(bayesplot)
library(bayesplot)
library(ggplot2)
library(posterior)

# After sampling:
# stan_out <- NCovid$sample(...)

# Extract posterior draws in array format
draws_array <- stan_out$draws(format = "array")

# Save permanently
#saveRDS(draws_array, "./output/Mhet_McMc_array.rds")
#str(readRDS("./output/Mhet_McMc_array.rds"))


# Load draws and diagnostics
draws_array <- readRDS("./output/Mhet_McMc_array.rds")

#remove NAS
param_names <- dimnames(draws_array)$variable

na_params <- sapply(param_names, function(p) {
  any(is.na(draws_array[, , p]))
})

clean_draws <- draws_array[, , !na_params]


params_to_check <- c(
  "beta[3]", "v[3]", "I0[3]", "p_reported",
  "beta[10]", "v[10]", "I0[10]", "p_reported"
)

# Traceplots
for (param in params_to_check) {
  trace_plot <- mcmc_trace(clean_draws, pars = param)
  safe_name <- gsub("\\[|\\]", "_", param)
  ggsave(paste0("Traceplot_", safe_name, ".png"), trace_plot, width = 8, height = 5, dpi = 300)
}

# Density plots
for (param in params_to_check) {
  dens_plot <- mcmc_dens_overlay(clean_draws, pars = param)
  safe_name <- gsub("\\[|\\]", "_", param)
  ggsave(paste0("Density_", safe_name, ".png"), dens_plot, width = 8, height = 5, dpi = 300)
}

for (param in params_to_check) {
  trace_plot <- mcmc_trace(clean_draws, pars = param)
  print(trace_plot)   # <-- shows plot in RStudio Viewer
  
  safe_name <- gsub("\\[|\\]", "_", param)
  ggsave(paste0("Traceplot_", safe_name, ".png"), trace_plot, width = 8, height = 5, dpi = 300)
}

for (param in params_to_check) {
  dens_plot <- mcmc_dens_overlay(clean_draws, pars = param)
  print(dens_plot)   # <-- shows plot
  
  safe_name <- gsub("\\[|\\]", "_", param)
  ggsave(paste0("Density_", safe_name, ".png"), dens_plot, width = 8, height = 5, dpi = 300)
}

print(dens_plot)
print(trace_plot)
ggsave("Output/Figure/Figure_dens_plot.png", width = 7, height = 5, bg = "white", dpi=1000)


#posterior pair plot
library(bayesplot)

edo_params <- c("beta[3]", "v[3]", "I0[3]", "p_reported")

pair_plot_edo <- mcmc_pairs(clean_draws, pars = edo_params)
print(pair_plot_edo)
ggsave("Pairs_Edo.png", pair_plot_edo, width = 10, height = 10, dpi = 300)


rivers_params <- c("beta[10]", "v[10]", "I0[10]", "p_reported")

pair_plot_rivers <- mcmc_pairs(clean_draws, pars = rivers_params)
print(pair_plot_rivers)
ggsave("Pairs_Rivers.png", pair_plot_rivers, width = 10, height = 10, dpi = 300)

##autocorrelatio plot
for (param in params_to_check) {
  ac_plot <- mcmc_acf(clean_draws, pars = param)
  print(ac_plot)
  
  safe_name <- gsub("\\[|\\]", "_", param)
  ggsave(paste0("ACF_", safe_name, ".png"), ac_plot, width = 8, height = 5, dpi = 300)
}


# pair plot is using posterior correlation others
edo_params <- c("beta[3]", "v[3]", "I0[3]", "p_reported")
rivers_params <- c("beta[10]", "v[10]", "I0[10]", "p_reported")
pair_plot_edo <- mcmc_pairs(clean_draws, pars = edo_params)
print(pair_plot_edo)   #pair plot for Edo
ggsave("Pairs_Edo.png", pair_plot_edo, width = 10, height = 10, dpi = 300)
pair_plot_rivers <- mcmc_pairs(clean_draws, pars = rivers_params)
print(pair_plot_rivers)



#Quantitative diagnostics#####
library(posterior)
library(bayesplot)

# Extract posterior draws
#draws <- stan_out$draws()

# Extract sampler diagnostics
diag <- stan_out$sampler_diagnostics()

# Save both
#saveRDS(draws, "./Output/draws_Mhet_NB.rds")
#saveRDS(diag,  "./Output/diag_Mhet_NB.rds")

# ---------------------------------------------------------
# Load draws and diagnostics
# ---------------------------------------------------------
draws <- readRDS("./Output/draws_Mhet_NB.rds")
diag  <- readRDS("./Output/diag_Mhet_NB.rds")

# ---------------------------------------------------------
# Rhat and ESS
# ---------------------------------------------------------
summ <- summarise_draws(draws)

max_rhat     <- max(summ$rhat,      na.rm = TRUE)
min_ess_bulk <- min(summ$ess_bulk,  na.rm = TRUE)
min_ess_tail <- min(summ$ess_tail,  na.rm = TRUE)

cat("Maximum Rhat:        ", max_rhat, "\n")
cat("Minimum bulk ESS:    ", min_ess_bulk, "\n")
cat("Minimum tail ESS:    ", min_ess_tail, "\n")

# ---------------------------------------------------------
# Divergent transitions (3D array)

var_div <- which(dimnames(diag)$variable == "divergent__")
divergent_vals <- diag[, , var_div]
n_divergent <- sum(divergent_vals)
cat("Divergent transitions:", n_divergent, "\n")

# ---------------------------------------------------------
# E-BFMI (manual computation)
# ---------------------------------------------------------
var_energy <- which(dimnames(diag)$variable == "energy__")
energy_vals <- diag[, , var_energy]

# Flatten all chains
energy_vec <- as.vector(energy_vals)

# BFMI formula: Var(E) / mean(diff(E)^2)
bfmi_val <- var(energy_vec) / mean(diff(energy_vec)^2)

cat("E-BFMI:               ", bfmi_val, "\n")

#end diagnostic

#posterior check and posterior mean curve######
##with legend

library(posterior)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(ggplot2)

#Extract posterior draws for pred_cases
draws_pred <- as_draws_df(stan_out$draws("pred_cases"))

n_states <- length(t_last)
max_days <- dim(cases_matrix)[1]

t_last_df <- tibble(state = 1:n_states, t_last = t_last)

#Convert pred_cases[s,i] into long format
long_pred <- draws_pred %>%
  pivot_longer(
    cols = starts_with("pred_cases["),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(
    var   = str_remove_all(var, "pred_cases\\[|\\]"),
    state = as.integer(str_extract(var, "^[0-9]+")),
    day   = as.integer(str_extract(var, "(?<=,)[0-9]+"))
  ) %>%
  left_join(t_last_df, by = "state") %>%
  filter(day <= t_last)

# 3. Posterior predictive mean + 95% intervals
summary_pred <- long_pred %>%
  group_by(state, day) %>%
  summarise(
    mean   = mean(value),
    median = median(value),
    lower  = quantile(value, 0.025),
    upper  = quantile(value, 0.975),
    .groups = "drop"
  )

# 4. Observed data
cases_subset <- as.matrix(cases_matrix[1:max(t_last), ])

obs_df <- map_dfr(1:n_states, function(s) {
  tibble(
    state    = s,
    day      = 1:(t_last[s] - 1),
    Observed = cases_subset[1:(t_last[s] - 1), s]
  )
})

#  Join predictions + observations
plot_df <- summary_pred %>%
  left_join(obs_df, by = c("state", "day"))

state_labels <- setNames(selected_states, 1:n_states)

# Build PPC plot WITH LEGEND
p_ppc <- ggplot(plot_df, aes(x = day)) +
  
  # 95% interval ribbon (legend via fill)
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = "95% Interval"),
    alpha = 0.35
  ) +
  
  # Posterior mean curve (legend via color)
  geom_line(
    aes(y = mean, color = "Posterior Mean"),
    size = 1
  ) +
  
  # Observed data (legend via color)
  geom_point(
    aes(y = Observed, color = "Observed Data"),
    size = 0.8
  ) +
  
  facet_wrap(~state, labeller = as_labeller(state_labels), scales = "free_y") +
  
  scale_fill_manual(
    name = "Posterior Predictive",
    values = c("95% Interval" = "lightblue")
  ) +
  
  scale_color_manual(
    name = "Posterior Predictive",
    values = c(
      "Posterior Mean" = "blue",
      "Observed Data"  = "black"
    )
  ) +
  
  labs(
    x = "Day",
    y = "Cases",
    title = "Posterior Predictive Intervals and Mean Curve by State"
  ) +
  
 # theme_minimal(base_size = 12) +
  theme_bw(base_size = 12)
  theme(legend.position = "bottom")

# 7. Show the plot
print(p_ppc)

# 8. Save the plot
ggsave("Posterior_Predictive_Check_legend.png", p_ppc, width = 12, height = 10, dpi = 600)

#posterior predict


#prior vs posterior used for evidence of identifiability######
#save in
# Extract posterior draws
#draws <- stan_out$draws()

# Save draws
#saveRDS(draws, "./Output/draws_model11_heterogeneity_NB.rds")

# Save metadata needed for plotting
#saveRDS(t_last, "./Output/t_last.rds")
#saveRDS(cases_matrix, "./Output/cases_matrix.rds")
#saveRDS(data_sir, "./Output/data_sir_Mhet_prior.rds")


rm(list = ls())

#load in
#stan_out <- readRDS("./Output/stan_out_Mhet_prior.rds")
#data_sir <- readRDS("./Output/data_sir_Mhet_prior.rds")
draws <- readRDS("./Output/draws_model11_heterogeneity_NB.rds")
t_last <- readRDS("./Output/t_last.rds")
cases_matrix <- readRDS("./Output/cases_matrix.rds")
data_sir <- readRDS("./Output/data_sir_Mhet_prior.rds")

#convert to draws
draws_df <- as_draws_df(draws)


v3_post  <- draws_df[["v[3]"]]
v10_post <- draws_df[["v[10]"]]

beta3_post  <- draws_df[["beta[3]"]]
beta10_post <- draws_df[["beta[10]"]]

p_reported_post <- draws_df[["p_reported"]]

I0_edo_post    <- draws_df[["I0[3]"]]
I0_rivers_post <- draws_df[["I0[10]"]]

R0_post <- draws_df[["R0"]]


#Prior vesrsus posterior#####
#Number of prior samples
n_prior <- 5000

# Priors for mixture parameters
theta_ss_prior <- rbeta(n_prior, 1, 1)
slab_sd_prior  <- rgamma(n_prior, 1, 10)

# spike_sd comes from your Stan data list
spike_sd <- data_sir$spike_sd

# Sample v from spike-and-slab prior
v_prior <- numeric(n_prior)
for (i in 1:n_prior) {
  if (runif(1) < theta_ss_prior[i]) {
    v_prior[i] <- rnorm(1, 0, spike_sd)      # spike
  } else {
    v_prior[i] <- rnorm(1, 0, slab_sd_prior[i])  # slab
  }
}


##combine prior plus posteriior
df <- data.frame(
  value = c(v_prior, v3_post, v10_post),
  type = c(
    rep("Prior (Spike-and-Slab)", length(v_prior)),
    rep("Posterior v[3] (Edo)", length(v3_post)),
    rep("Posterior v[10] (Rivers)", length(v10_post))
  )
)


#plot


p <- ggplot(df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.45) +
  labs(
    title = "Prior vs Posterior for v Parameter",
    x = "v value",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

print(p)

ggsave("prior_posterior_v.png", p, width = 8, height = 5, dpi = 300)

#for beta
library(posterior)
library(ggplot2)

#extract draws
#draws <- stan_out$draws()
#saveRDS(draws, "./Output/draws_model11_heterogeneity_NB.rds")

#read in
# Load saved posterior draws
draws <- readRDS("./Output/draws_model11_heterogeneity_NB.rds")

# Convert to draws_df
draws_df <- as_draws_df(draws)

# 1. Extract posterior draws

beta3_post  <- as_draws_df(stan_out$draws("beta[3]"))$`beta[3]`
beta10_post <- as_draws_df(stan_out$draws("beta[10]"))$`beta[10]`

# 2. Simulate prior samples
# Stan prior: beta ~ normal(0.5, 0.1)

n_prior <- 5000
beta_prior <- rnorm(n_prior, mean = 0.5, sd = 0.1)


# 3. Combine prior  and posterior

df_beta <- data.frame(
  value = c(beta_prior, beta3_post, beta10_post),
  type = c(
    rep("Prior (Normal 0.5, 0.1)", length(beta_prior)),
    rep("Posterior beta[3] (Edo)", length(beta3_post)),
    rep("Posterior beta[10] (Rivers)", length(beta10_post))
  )
)

p_beta <- ggplot(df_beta, aes(x = value, fill = type)) +
  geom_density(alpha = 0.45) +
  labs(
    title = "Prior vs Posterior for beta Parameter",
    x = "beta value",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

print(p_beta)

ggsave("prior_posterior_beta.png", p_beta, width = 8, height = 5, dpi = 300)


#for all
#Prior vs posterior fors######
library(posterior)
library(ggplot2)

# Choose two states to illustrate (Edo = 3, Rivers = 10)
state_edo    <- 3
state_rivers <
  # Extract posterior draws
  
  draws_df <- as_draws_df(stan_out$draws())

post <- list(
  beta_edo    = draws_df[[sprintf("beta[%d]", state_edo)]],
  beta_rivers = draws_df[[sprintf("beta[%d]", state_rivers)]],
  v_edo       = draws_df[[sprintf("v[%d]", state_edo)]],
  v_rivers    = draws_df[[sprintf("v[%d]", state_rivers)]],
  p_reported  = draws_df[["p_reported"]],
  I0_edo      = draws_df[[sprintf("I0[%d]", state_edo)]],
  I0_rivers   = draws_df[[sprintf("I0[%d]", state_rivers)]],
  R0          = draws_df[["R0"]]
)


# Simulate priors

n_prior <- 5000

prior <- list(
  beta        = rnorm(n_prior, 0.5, 0.1),
  p_reported  = rbeta(n_prior, 0.06, 8)
)

# spike-and-slab prior for v
theta_ss_prior <- rbeta(n_prior, 1, 1)
slab_sd_prior  <- rgamma(n_prior, 1, 10)
spike_sd       <- data_sir$spike_sd

v_prior <- numeric(n_prior)
for (i in 1:n_prior) {
  if (runif(1) < theta_ss_prior[i]) {
    v_prior[i] <- rnorm(1, 0, spike_sd)
  } else {
    v_prior[i] <- rnorm(1, 0, slab_sd_prior[i])
  }
}
prior$v <- v_prior

# I0 prior depends on p_reported prior
prior$I0 <- rnorm(n_prior, mean = 1 / prior$p_reported,
                  sd   = 1 / prior$p_reported)

# R0 prior from beta prior
prior$R0 <- prior$beta / 0.2

# -----------------------------
#plot n save
plot_prior_posterior <- function(prior_vec, post_edo_vec, post_rivers_vec,
                                 label) {
  
  df <- data.frame(
    value = c(prior_vec, post_edo_vec, post_rivers_vec),
    type  = c(
      rep("Prior", length(prior_vec)),
      rep(sprintf("Posterior %s (Edo)", label),    length(post_edo_vec)),
      rep(sprintf("Posterior %s (Rivers)", label), length(post_rivers_vec))
    )
  )
  
  p <- ggplot(df, aes(x = value, fill = type)) +
    geom_density(alpha = 0.45) +
    labs(
      title = sprintf("Prior vs Posterior for %s", label),
      x = sprintf("%s value", label),
      y = "Density"
    ) +
    theme_minimal(base_size = 14)
  
  print(p)  # show in RStudio
  
  ggsave(
    filename = paste0("prior_posterior_", label, ".png"),
    plot = p,
    width = 8, height = 5, dpi = 300
  )
}


#Generate and save all plots

plot_prior_posterior(prior$beta, post$beta_edo, post$beta_rivers, "beta")
plot_prior_posterior(prior$v,    post$v_edo,    post$v_rivers,    "v")
plot_prior_posterior(prior$p_reported, post$p_reported, post$p_reported, "p_reported")
plot_prior_posterior(prior$I0,  post$I0_edo,   post$I0_rivers,    "I0")
plot_prior_posterior(prior$R0,  post$R0,       post$R0,           "R0")

