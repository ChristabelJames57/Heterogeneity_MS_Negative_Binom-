#Abstract and disscusion computation#########

#HIT computation FOR MHOM AND MHET in the abstract#####

library(dplyr)

# Fixed parameter
gamma <- 0.2

# Inputs
beta_vec <- c(#mhet beta
  0.3009, 0.2935, 0.2676, 0.2896, 0.2367,
  0.2656, 0.2914, 0.2689, 0.2565, 0.2506
)

v_vec <- c(#Mhet v
  3.0728, 2.5497, 1.7406, 3.1850, 1.8934,
  2.2834, 2.8745, 3.3145, 1.7874, 2.8323
)

# state data frame
state_df <- data.frame(
  state = factor(1:10),
  R0    = beta_vec / gamma,
  v_hat = v_vec
)

#COMPUTE HIT#####

hit_summary <- state_df %>%
  mutate(
    HIT_hom = 1 - 1 / R0,
    HIT_het = 1 - R0^(-1 / (1 + v_hat^2)),
    percent_drop = 100 * (HIT_hom - HIT_het) / HIT_hom
  )

print(hit_summary)

#mean HIT
mean_HIT_hom <- mean(hit_summary$HIT_hom)
mean_HIT_het <- mean(hit_summary$HIT_het)

#% mean HIT

mean_HIT_hom_percent <- 100 * mean_HIT_hom
mean_HIT_het_percent <- 100 * mean_HIT_het

# Mean %drop
mean(hit_summary$percent_drop)
#end


#Mean attack rate and Final epidemic size#####
library(deSolve)
library(tidyverse)
library(purrr)

#### PARAMETERS
gamma <- 0.2

pop_order <- c(
  Delta  = 5663362,
  Ebonyi = 2880383,
  Edo    = 4235595,
  Enugu  = 4411119,
  Kaduna = 8252366,
  Kwara  = 3192893,
  Ondo   = 4671695,
  Osun   = 4705589,
  Oyo    = 7840864,
  Rivers = 7303924
)

pN <- c(
  Delta  = 0.3756188,
  Ebonyi = 0.4886132,
  Edo    = 0.8194811,
  Enugu  = 0.3356231,
  Kaduna = 0.6693090,
  Kwara  = 0.5351148,
  Ondo   = 0.4057060,
  Osun   = 0.2896763,
  Oyo    = 0.7841961,
  Rivers = 0.3727742
)

times <- seq(0, 600, by = 1)

beta_vec <- c(
  0.3009, 0.2935, 0.2676, 0.2896, 0.2367,
  0.2656, 0.2914, 0.2689, 0.2565, 0.2506
)

v_vec <- c(
  3.0728, 2.5497, 1.7406, 3.1850, 1.8934,
  2.2834, 2.8745, 3.3145, 1.7874, 2.8323
)

beta_vec_Mhom_het <- c(
  0.29932, 0.29263, 0.2679, 0.2820, 0.2367,
  0.26533, 0.29025, 0.2681, 0.2567, 0.25006
)

v_vec_Mhom_het <- c(
  2.1308, 2.0596, 1.97007, 2.0675, 1.8652,
  1.9491, 2.0815, 1.9782, 1.9611, 1.9512
)

I0_vec_Mhet <- round(c(
  2176.529, 3208.261, 2374.598, 3019.818, 5257.560,
  2121.510, 1229.894, 2663.764, 6255.512, 9235.500
))

I0_vec_Mhom <- round(c(
  320.851, 659.3278, 494.2246, 584.7005, 1266.619,
  489.3966, 194.1634, 595.5991, 1337.0867, 1341.4041
))

I0_vec_Mhom_het <- round(c(
  447.1475, 949.2842, 732.2879, 831.2845, 1855.078,
  1124.465, 272.8941, 858.3466, 1971.380, 1944.409
))

#### MODELS
SIR_hetero <- function(time, state, parameters) {
  S <- state[1]; I <- state[2]; R <- state[3]
  beta  <- parameters["beta"]
  gamma <- parameters["gamma"]
  N     <- parameters["N"]
  v     <- parameters["v"]
  frac <- max(min(S / N, 1), 1e-8)
  dS <- -beta * I * frac^(1 + v^2)
  dI <-  beta * I * frac^(1 + v^2) - gamma * I
  dR <-  gamma * I
  list(c(dS, dI, dR))
}

SIR_homo <- function(time, state, parameters) {
  S <- state[1]; I <- state[2]; R <- state[3]
  beta  <- parameters["beta"]
  gamma <- parameters["gamma"]
  N     <- parameters["N"]
  frac <- max(min(S / N, 1), 1e-8)
  dS <- -beta * I * frac
  dI <-  beta * I * frac - gamma * I
  dR <-  gamma * I
  list(c(dS, dI, dR))
}

simulate_final_size <- function(N, model, beta, gamma, I0, v = NULL) {
  I0 <- min(I0, 0.99 * N)
  init <- c(S = N - I0, I = I0, R = 0)
  parms <- c(beta = beta, gamma = gamma, N = N)
  if (!is.null(v)) parms <- c(parms, v = v)
  
  out <- ode(y = init, times = times, func = model, parms = parms)
  out_df <- as.data.frame(out)
  
  S_T <- tail(out_df$S, 1)
  final_size <- N - S_T
  
  tibble(final_size = final_size)
}

#### RUN SIMULATION
results <- map_dfr(seq_along(pop_order), function(i) {
  
  state <- names(pop_order)[i]
  N <- pop_order[[i]]
  pN_i <- pN[[state]]
  
  tibble(
    final_size_Mhet = simulate_final_size(N, SIR_hetero, beta_vec[i], gamma, I0_vec_Mhet[i], v_vec[i])$final_size,
    final_size_Mhom = simulate_final_size(N, SIR_homo, beta_vec[i], gamma, I0_vec_Mhom[i])$final_size,
    final_size_Mhom_pi = simulate_final_size(N * pN_i, SIR_homo, beta_vec[i], gamma, I0_vec_Mhom[i])$final_size,
    final_size_Mhom_het = simulate_final_size(N * pN_i, SIR_hetero, beta_vec_Mhom_het[i], gamma, I0_vec_Mhom_het[i], v_vec_Mhom_het[i])$final_size
  )
})

N_mean <- mean(pop_order)

attack_rate_table_pretty <- tibble(
  model = c("Mhet", "Mhom", "Mhom_pi", "Mhom_het"),
  
  mean_final_size = c(
    mean(results$final_size_Mhet),
    mean(results$final_size_Mhom),
    mean(results$final_size_Mhom_pi),
    mean(results$final_size_Mhom_het)
  ),
  
  mean_attack_rate = 100 * c(
    mean(results$final_size_Mhet) / N_mean,
    mean(results$final_size_Mhom) / N_mean,
    mean(results$final_size_Mhom_pi) / N_mean,
    mean(results$final_size_Mhom_het) / N_mean
  ),
  
  range = c(
    paste0(round(min(results$final_size_Mhet / N_mean)*100,2), "–", round(max(results$final_size_Mhet / N_mean)*100,2), "%"),
    paste0(round(min(results$final_size_Mhom / N_mean)*100,2), "–", round(max(results$final_size_Mhom / N_mean)*100,2), "%"),
    paste0(round(min(results$final_size_Mhom_pi / N_mean)*100,2), "–", round(max(results$final_size_Mhom_pi / N_mean)*100,2), "%"),
    paste0(round(min(results$final_size_Mhom_het / N_mean)*100,2), "–", round(max(results$final_size_Mhom_het / N_mean)*100,2), "%")
  )
)

attack_rate_table_pretty
#end attack rate and final epidemic size 


#Correlation of all parameters#####
library(dplyr)
library(tidyr)
library(stringr)
library(posterior)

# Extract all posterior draws into a tidy data frame
posterior_df <- as_draws_df(stan_out$draws())

# Save to RDS
saveRDS(posterior_df, "./output/posterior_all_main_parameters2.rds")

#read in
draws <- readRDS("./output/posterior_all_main_parameters2.rds")

# Convert to long format (NOW INCLUDING p_reported)
draws_long <- draws %>%
  pivot_longer(
    cols = matches("^(v|p_reported|beta|R0|I0)\\["),
    names_to = c("param", "state"),
    names_pattern = "(v|p_reported|beta|R0|I0)\\[(\\d+)\\]",
    values_to = "value"
  ) %>%
  mutate(state = as.integer(state))


### RESHAPE DRAWS PER STATE
draws_wide <- draws_long %>%
  pivot_wider(
    names_from = param,
    values_from = value
  )


## COMPUTE CORRELATION PER STATE
state_correlations <- draws_wide %>%
  group_by(state) %>%
  summarise(
    cor_v_I0        = cor(v, I0),
    cor_v_beta      = cor(v, beta),
    cor_v_R0        = cor(v, R0),
    cor_v_prep      = cor(v, p_reported),
    
    cor_I0_beta     = cor(I0, beta),
    cor_I0_R0       = cor(I0, R0),
    cor_I0_prep     = cor(I0, p_reported),
    
    cor_beta_R0     = cor(beta, R0),
    cor_beta_prep   = cor(beta, p_reported),
    cor_R0_prep     = cor(R0, p_reported)
  )

print(state_correlations)


## IDENTIFY STRONG CORRELATION
strong_corr <- state_correlations %>%
  pivot_longer(-state, names_to = "pair", values_to = "correlation") %>%
  filter(abs(correlation) > 0.6)

print(strong_corr)
print(strong_corr, n = Inf)#more rows
#end

