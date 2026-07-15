#THE Mhom-het model#####
library(bayesplot)
library(cmdstanr)
library(tidybayes)
library(tidyverse)
library(gridExtra)
library(outbreaks)
library(readr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(loo)
library(posterior)
library(writexl)
library(stringr)
library(purrr)


rm(list = ls())
set.seed(1275)

MAX_TIME <- 189

#selected states
selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu","Kaduna", "Kwara", "Ondo", "Osun","Oyo", "Rivers")

n_states <- length(selected_states)

pN <- c(
  0.3756188,
  0.4886132,
  0.8194811,
  0.3356231,
  0.6693090,
  0.5351148,
  0.4057060,
  0.2896763,
  0.7841961,
  0.3727742
)

##read in
df_N <- read.csv("./Data/Nigeria17.csv")

populations <- c(
  Ondo = 4671695, Edo = 4235595, Osun = 4705589,
  Ogun = 5217716, Oyo = 7840864, Lagos = 12550598,
  Kano = 13076892, Kwara = 3192893, Delta = 5663362,
  Ebonyi = 2880383, Enugu = 4411119, Kaduna = 8252366,
  Bayelsa = 2277961, Benue = 5741815, Fct = 3564126,
  Gombe = 3256962, Rivers = 7303924, Bauchi = 6537314,
  `Cross River` = 3866269
)

end_days <- c(
  Delta = 140, Ebonyi = 119, Edo = 189, Enugu = 128,
  Kaduna = 200, Kwara = 154, Ondo = 152, Osun = 127,
  Oyo = 174, Rivers = 157
)

t_last <- end_days[selected_states]

#CASE MATRIX
cases_matrix <- df_N %>%
  select(all_of(selected_states)) %>%
  as.matrix()

##COMPUTE EFFECTIVE POPUL
N_raw <- populations[selected_states]

N_effective <- as.integer(round(pN * N_raw))

print(data.frame(
  State = selected_states,
  N_raw = N_raw,
  pN = pN,
  N_effective = N_effective
))

data_sir <- list(
  max_days = max(t_last),
  n_states = n_states,
  t0 = 0,
  t_last = t_last,
  N = N_effective,      # N = pN × N
  cases = cases_matrix[1:max(t_last), ],
  spike_sd = 1e-3
)

#compile
NCovid <- cmdstan_model("./Scripts/Mhom_het_model_NB.stan")

iters <- 2000

stan_out <- NCovid$sample(data = data_sir,iter_warmup = iters, iter_sampling = iters,parallel_chains = 4, seed = 123,adapt_delta = 0.95)

print(stan_out, variables = c( "MeanV", "meanbeta", "meanp_reported", "meanI0", "meanR0",  "theta_ss", "v", "beta", "I0", "phi" ),digits = 6, max_rows = 60)


#----Model selction----
log_lik_array <- stan_out$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik_array)
waic_result <- waic(log_lik_array)

print(loo_result)
print(waic_result)

stan_out$diagnostic_summary()#Check divergence


# Extract posterior draws for prediction plot
draws_df <- as_draws_df(stan_out$draws("pred_cases"))

# Convert to long format
long_draws <- draws_df %>%
  pivot_longer(
    cols = starts_with("pred_cases["),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(
    var = str_remove_all(var, "pred_cases\\[|\\]"),
    state = as.integer(str_extract(var, "^[0-9]+")),
    day = as.integer(str_extract(var, "(?<=,)[0-9]+"))
  )

# Keep only valid days per state
t_last_df <- tibble(state = 1:length(t_last), t_last = t_last)

long_draws <- long_draws %>%
  left_join(t_last_df, by = "state") %>%
  filter(day < t_last)

# Posterior median only (NO CI)
summary_df <- long_draws %>%
  group_by(state, day) %>%
  summarise(
    median = median(value),
    .groups = "drop"
  )

# Observed data
cases_subset <- cases_matrix[1:max(t_last), ]

obs_df <- map_dfr(1:length(t_last), function(s) {
  tibble(
    state = s,
    day = 1:(t_last[s] - 1),
    Observed = cases_subset[1:(t_last[s] - 1), s]
  )
})

# Merge
plot_df <- summary_df %>%
  left_join(obs_df, by = c("state", "day"))

state_labels <- setNames(selected_states, 1:length(selected_states))

#plot
ggplot(plot_df, aes(x = day)) +
  geom_point(aes(y = Observed),
             color = "black",
             size = 0.8) +
  geom_line(aes(y = median),
            color = "blue",
            size = 1) +
  facet_wrap(~state,
             labeller = as_labeller(state_labels),
             scales = "free_y") +
  labs(
    x = "Day",
    y = "Cases",
   # title = "Observed vs Predicted Cases by State"
  ) +
  theme_minimal(base_size = 12)

###TO new SAVE MODEL CODE ESTIMATING N=N*PN to use for prediction
draws_df_m44 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m44, file = "draws_model44_PNheterogeneity_NB.rds")
#end