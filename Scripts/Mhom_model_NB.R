#Mhom model using Negative Binom observation error#####
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
library(deSolve)
rm(list = ls()) #clears the memory

##FIXING V (OR NO V) AND ESTIMATIN PN FOR REAL ((S / (pN*N))) NO DECAY paramter W

set.seed(1275) 
total_number_of_time_series=19#for simulation
MAX_TIME<-189#IF real
#MAX_TIME<-100# if trying simulation
N=matrix(data=0, nrow = total_number_of_time_series)
tru_incidence_matrix <- obs_incidence_matrix<- matrix(data=0, nrow = MAX_TIME , ncol = total_number_of_time_series)#all in a matrix with that row n column
Incidence=matrix(data=0, nrow = MAX_TIME,ncol=total_number_of_time_series)#so we have day row n outbreak col
obs_incidence_matrix2 <- matrix(0, nrow = MAX_TIME, ncol = total_number_of_time_series)




##Using ((S / (pN*N)))  Here pN=1 if simualtion bc not estimating pN. if not simualtion, use pN here to estimate
SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma  # recompute beta from R0 and gamma
    dS <- -beta * I * (S /(pN* N[i]))# bc estimating PN for real and not for simula
    dI <-  beta * I * (S / (pN*N[i])) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}


cmdstanr_max_rows=15
max_rows = 15
#setwd("~/Desktop/Dans Folder 2023/People/Christabel James/Stan code/Attempt3")

The_mode="R"#change only here from  S to R N chge prior for rp
if (The_mode=="R"){
  df_N <- read.csv("./Data/Nigeria17.csv") # Columns: Day, Edo,Lagos, Ondo, Ogun, Osun, Oyo, 
  populations <- c(Ondo = 4671695, Edo = 4235595, Osun = 4705589, Ogun = 5217716, Oyo = 7840864, Lagos = 12550598, Kano = 13076892, Kwara = 3192893, Delta = 5663362, Ebonyi = 2880383, Enugu = 4411119, Kaduna = 8252366, Bayelsa = 2277961, Benue = 5741815, Fct = 3564126, Gombe = 3256962, Rivers = 7303924, Bauchi = 6537314, `Cross River` = 3866269)
  #populations <- c(Ondo = 4671695, Edo = 4235595, Osun = 4705589, Ogun = 5217716, Oyo = 7840864, Lagos = 12550598, Kano = 13076892, Kwara = 3192893, Delta = 5663362, Ebonyi = 2880383, Enugu = 4411119, Kaduna = 8252366)
  #end_days<-c(Ondo = 152, Edo = 189, Osun = 127, Ogun = 197, Oyo = 174, Lagos = 228,Kano =144, Kwara = 154, Delta = 140, Ebonyi = 119, Enugu = 128, Kaduna = 200)#real data 10 states
  end_days<-c(Ondo = 152, Edo = 189, Osun = 127, Ogun = 197, Oyo = 174, Lagos = 228,Kano =144, Kwara = 154, Delta = 140, Ebonyi = 119, Enugu = 128, Kaduna = 200, Bayelsa =89, Benue =86, Fct= 225, Gombe = 143, Rivers= 157, Bauchi = 129, `Cross River` =89)#real data 15 states
  
  
  # Select States (ten)
  selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara","Ondo", "Osun", "Oyo", "Rivers")#arrge alphabe 13 no Benue
  
  # Prepare Data 
  n_states <- length(selected_states)
  #n_days <- length(df_N$Day)
  t_last <- end_days[selected_states]#end days for selcted states alone
  
  cases_matrix <- df_N %>%#matrix of case data
    select(all_of(selected_states)) %>%#selects the column for specified states
    as.matrix(nrow=t_last,ncol=n_states)#converts to matrix
  
  
  N_vector <- populations[selected_states]
  #end real data
  
} else{#begin simulated data
  
  #beta <- 0.3765#Edo real data
  beta <- 0.4#
  gamma <- 0.2
  #v <- 1.4946 # Silence v in true bc not there
  #pN<-0.4#ONLY SILENCE IF USING SIMULATIO. IF REAL, NO TRUTH 
  R0 <- beta / gamma
  #I0 <- 5
  reporting_prob <- 0.0075#using now
  I0 <- rpois(1,1/ reporting_prob)#DIFFERENT Initial conditions
  
  sim_r <- 2.1
  noof_outbreaks<-5#as you want but must be less than population 12
  n_states <- noof_outbreaks
  selected_states <- paste0("SimState_", 1:noof_outbreaks)
  N_vector<-rep(0, each = noof_outbreaks)
  t_last<-rep(0,each = noof_outbreaks)
  #N_vector<-matrix(0,nrow=noof_outbreaks)# for simu
  #t_last<-matrix(0, nrow = noof_outbreaks)
  
  
  for (i in 1:noof_outbreaks) {
    N[i] <- 4235595#Edo
    S0 <- N[i] - I0
    R0_init <- 0
    y_init <- c(S = S0, I = I0, R = R0_init)
    
    times <- seq(0, MAX_TIME, by = 1)
    
    
    parameters <- c(beta = beta, gamma = gamma, pN=pN)#remove pN if not estimating it that is for simualtion. Use her if Real
    
    #parameters <- c(beta = beta, gamma = gamma)# Use like this if simulation
    
    out <- ode(y = y_init, times = times, func = SIRode, parms = parameters)
    out_df <- as.data.frame(out)
    
    
    incidence <- -diff(out_df$S)# Calculate incidence as new infections (delta S)
    incidence[incidence < 0] <- 0 
    
    #True process for incidence (deltaS)
    tru_incidence_matrix[,i ] <- incidence
    
    #Observed process
    obs_incidence_matrix[i, ]<-rnbinom(MAX_TIME, mu=reporting_prob* Incidence[i,],size=sim_r)#-VE BINOM
   # obs_incidence_matrix[,i ] <- rnbinom(MAX_TIME, mu = reporting_prob * incidence + 1e-6, size = sim_r)#incidence computd from ODe
    obs_incidence_matrix2[, i ] <- obs_incidence_matrix[,i ]#so we have day row n outbreak col
    N_vector[i] <- N[i]
    t_last[i]<-MAX_TIME
    
  }
  
  cases_matrix <- obs_incidence_matrix[,1:noof_outbreaks] #matrix of 4 cols
  
}#end simulated

#I0 <- 1
#y0_array <- t(sapply(N_vector, function(N) c(N - I0, I0, 0)))  # [n_states, 3]#here R=0, I=1 and S=N-1

data_sir <- list(max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,N = N_vector,cases = cases_matrix[1:max(t_last), ])#pas as input to stan NO PHI N NO IO HERE BUT STAN

#NCovid <-cmdstan_model("Mhom_model_NB.stan")#NEW ESTI_Feb_manuscrip same i0 prior
NCovid <-cmdstan_model("./Scripts/Mhom_model_NB.stan")#

iters=2000

##Regularises where the initial value begins
init_fun <- function() {
  list(
    pN = rep(0.6, 10),               # vector of same length
    beta = rep(0.4, 10),       # now a vector of length n_states
    p_reported = 0.02,                   # scalar is fine (global parameter)
    I0 = rep(1/0.02, 10)  # vector of same length
  )
}

stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4,init = init_fun)
#stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4, seed = 0)

print(stan_out, variables = c("MeanpN","meanbeta", "meanR0", "meanp_reported", "meanI0", "beta", "R0","pN", "I0", "phi"), digits = 8,max_rows = 70)#po#summa if simualtion, Now I0


##.... to check for R0 CI##
draws <- stan_out$draws()
R0_samples <- as_draws_df(draws)[["R0"]]
quantile(R0_samples, probs = c(0.025, 0.5, 0.975))#2.5% and 7.5% and mean 50%

#----Model selction----
log_lik_array <- stan_out$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik_array)
waic_result <- waic(log_lik_array)

print(loo_result)
print(waic_result)

stan_out$diagnostic_summary()#Check divergence

#posterior_samples <- as.data.frame(stan_out$draws(format = "draws_df"))
posterior_samples <- stan_out$draws(format = "draws_df")#trace plots
#######################. code for traces ######################
all_vars <- variables(stan_out$draws())

# Exclude generated quantities with many entries or possible NAs
exclude_patterns <- c("pred_cases", "log_lik", "y", "incidence")
valid_vars <- all_vars[!sapply(all_vars, function(v) any(sapply(exclude_patterns, grepl, v)))]

# Split into chunks of size 10
chunk_size <- 2
var_chunks <- split(valid_vars, ceiling(seq_along(valid_vars) / chunk_size))

# Plot each chunk
for (i in seq_along(var_chunks)) {
  message(paste("Plotting trace for chunk", i, "/", length(var_chunks)))
  
  current_chunk <- var_chunks[[i]]
  
  # Try plotting and catch any errors due to NAs or missing variables
  tryCatch({
    p <- mcmc_trace(stan_out$draws(current_chunk), 
                    facet_args = list(ncol = 2)) +
      ggtitle(paste("Trace Plots: Chunk", i))
    print(p)
  }, error = function(e) {
    message("Error in chunk ", i, ": ", conditionMessage(e))
  })
}


#  Posterior vs Prior (not fully working)
pars <- c("beta", "pN", "phi_inv", "p_reported")
n <- 10000
prior <- tibble(
  beta = abs(rnorm(n * n_states, 0.5, 0.1)),
  pN = rbeta(n * n_states, 4, 4),
  phi= rexp(n * n_states, 50),
  #p_reported = rbeta(n * n_states, 1, 10000)#real data
  p_reported = rbeta(n * n_states, 1, 1000)#simulated data
)

prior_df <- prior %>%
  mutate(state = rep(selected_states, each = n)) %>%
  pivot_longer(cols = -state, names_to = "parameter", values_to = "value") %>%
  mutate(type = "Prior")

posterior_df <- posterior_samples %>%
  select(matches(paste0("^(", paste(pars, collapse = "|"), ")\\["))) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  separate(param, into = c("parameter", "index"), sep = "\\[", remove = TRUE) %>%
  mutate(index = str_remove(index, "\\]")) %>%
  mutate(state = selected_states[as.numeric(index)], type = "Posterior")

combined_df <- bind_rows(prior_df, posterior_df) %>%
  mutate(parameter = factor(parameter, levels = pars))

##pLot of posterior vs prior
ggplot(combined_df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_grid(parameter ~ state, scales = "free") +
  scale_fill_manual(values = c("Prior" = "gray70", "Posterior" = "blue")) +
  theme_minimal() +
  labs(title = "Prior vs Posterior Distributions", x = "Value", y = "Density", fill = NULL)



#######################. code for predictions ######################
# Extract posterior draws
draws_df <- as_draws_df(stan_out$draws("pred_cases"))

# Dimensions
n_states <- length(t_last)
max_days <- dim(cases_matrix)[1]

# Name t_last for joining
t_last_df <- tibble(state = 1:n_states, t_last = t_last)

# Pivot draws to long format and parse indices
long_draws <- draws_df %>%
  pivot_longer(cols = starts_with("pred_cases["), names_to = "var", values_to = "value") %>%
  mutate(
    var = str_remove_all(var, "pred_cases\\[|\\]"),
    state = as.integer(str_extract(var, "^[0-9]+")),
    day = as.integer(str_extract(var, "(?<=,)[0-9]+"))
  ) %>%
  left_join(t_last_df, by = "state") %>%
  filter(day <= t_last[state])  # Keep only valid days per state

# Summarise draws: median and 95% credible interval
summary_df <- long_draws %>%
  group_by(state, day) %>%
  summarise(
    median = median(value),
    #lower = quantile(value, 0.25),   # 25% → for 50% CI
    #upper = quantile(value, 0.75),   # 75% → for 50% CI
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )

# Observed data prep
cases_subset <- as.matrix(cases_matrix[1:max(t_last), ])

obs_df <- purrr::map_dfr(1:n_states, function(s) {
  tibble(
    state = s,
    day = 1:(t_last[s] - 1),
    Observed = cases_subset[1:(t_last[s] - 1), s]
  )
})

# Join predictions and observations
plot_df <- summary_df %>%
  left_join(obs_df, by = c("state", "day"))
#state_labels <- setNames(selected_states, c(1, 2, 3, 4))#replac 1234 with abc for four states
state_labels <- setNames(selected_states, c(1, 2, 3, 4, 5,6,7,8,9,10))#replace 1234 with names if 10 states

# Plot using facets
ggplot(plot_df, aes(x = day)) +
  geom_point(aes(y = Observed), color = "black", size = 0.8) +
  geom_line(aes(y = median), color = "blue", size = 1) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.15)+#to lighten the colour for bette visbility
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +##CI
  #facet_wrap(~state, scales = "free_y") +
  facet_wrap(~state, labeller = as_labeller(state_labels), scales = "free_y") +#FORlabeling states
  labs(x = "Day", y = "Cases", title = "Observed vs Predicted Cases by State") +
  theme_minimal(base_size = 12)


###TO SAVE MODEL CODE


###TO new SAVE MODEL CODE ESTIMATING PN SAME prior
draws_df_m33 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m33, file = "draws_model33_noheterogeneity_NB.rds")
#end

###save t_last
t_last <- c(152, 189, 127, 197, 174, 228, 144, 154, 140, 119)

saveRDS(t_last, "t_last.rds")

##save case matrix only once
library(tidyverse)

df_N <- read.csv("Nigeria17.csv")

selected_states <- c(
  "Delta", "Ebonyi", "Edo", "Enugu", "Kaduna",
  "Kwara", "Ondo", "Osun", "Oyo", "Rivers"
)

cases_matrix <- df_N %>%
  select(all_of(selected_states)) %>%
  as.matrix()

saveRDS(cases_matrix, "cases_matrix.rds")
saveRDS(selected_states, "selected_states.rds")




#####----FOR FINAL EPIDEMIC SIZE##
###new
##### FINAL EPIDEMIC SIZE
final_epidemic_size_draws <- stan_out$draws("final_epidemic_size")
library(posterior)
final_epidemic_size_df <- as_draws_df(final_epidemic_size_draws)
apply(final_epidemic_size_df, 2, mean)
apply(final_epidemic_size_df, 2, quantile, probs = c(0.025, 0.5, 0.975))

##mean
state_names <- selected_states
mean_final_epidemic_size <- setNames(
  apply(final_epidemic_size_df, 2, mean),
  state_names
)

mean_final_epidemic_size

####mean posterior
overall_final_epidemic_size <- rowMeans(final_epidemic_size_df)

mean(overall_final_epidemic_size )# 146735
quantile(overall_final_epidemic_size, probs = c(0.025, 0.5, 0.975))



#end
library(posterior)

draws <- as_draws_df(stan_out)

final_size <- draws %>%
  select(starts_with("final_epidemic_size"))
#
final_size_summary <- final_size %>%
  pivot_longer(everything(),
               names_to = "state",
               values_to = "final_size") %>%
  group_by(state) %>%
  summarise(
    mean = mean(final_size),
    median = median(final_size),
    lower = quantile(final_size, 0.025),
    upper = quantile(final_size, 0.975)
  )

final_size_summary#end

##mean final
mean_final_size_draws <- final_size %>%
  mutate(mean_across_states = rowMeans(across(everything())))

mean_across_states_summary <- mean_final_size_draws %>%
  summarise(
    mean = mean(mean_across_states),
    median = median(mean_across_states),
    lower = quantile(mean_across_states, 0.025),
    upper = quantile(mean_across_states, 0.975)
  )

mean_across_states_summary


###attack rate
final_size_summary <- final_size_summary %>%
  mutate(state = selected_states)

library(posterior)

pN_draws <- as_draws_df(stan_out) %>%
  select(starts_with("pN["))

pN_mean <- colMeans(pN_draws)   # length = 10

N_df <- tibble(
  state = selected_states,
  N_eff = as.numeric(pN_mean) * data_sir$N
)

final_size_summary <- final_size_summary %>%
  mutate(
    N_eff = as.numeric(pN_mean) * data_sir$N,
    attack_rate = mean / N_eff
  )
final_size_summary

###mean attack rate
mean_attack_rate <- final_size_summary %>%
  summarise(mean_attack_rate = mean(attack_rate))

mean_attack_rate


######for 95Ci
library(posterior)
library(dplyr)

draws <- as_draws_df(stan_out)

# Final epidemic size draws (one column per state)
final_draws <- draws %>% 
  select(starts_with("final_epidemic_size["))

# pN draws (one column per state)
pN_draws <- draws %>% 
  select(starts_with("pN["))

##comput
# Turn N into a vector aligned with columns
N_vec <- data_sir$N

attack_rate_draws <- final_draws

for (s in seq_along(N_vec)) {
  attack_rate_draws[[s]] <-
    final_draws[[s]] / (pN_draws[[s]] * N_vec[s])
}

##
attack_rate_draws <- attack_rate_draws %>%
  mutate(mean_attack_rate = rowMeans(across(everything())))
###ci
mean_attack_rate_ci <- attack_rate_draws %>%
  summarise(
    mean   = mean(mean_attack_rate),
    median = median(mean_attack_rate),
    lower  = quantile(mean_attack_rate, 0.025),
    upper  = quantile(mean_attack_rate, 0.975)
  )

mean_attack_rate_ci
#end

##use for real data model 2 to estimate Pn
SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma  # recompute beta from R0 and gamma
    dS <- -beta * I * (S / (pN*N[i]))^(1 + v^2)
    dI <-  beta * I * (S / (pN*N[i]))^(1 + v^2) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}


# Load package model 1dder
library(deSolve)

# Parameters
N <- 4000000         # total population
gamma <- 0.2         # recovery rate
beta <- 0.41         # transmission rate
t <- seq(0, 200, by = 1)   # time steps

# Define the SIR model
SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS.dt <- -beta * (S/N)^p * I
    dI.dt <-  beta * (S/N)^p * I - gamma * I
    dR.dt <-  gamma * I
    return(list(c(dS.dt, dI.dt, dR.dt)))
  })
}

# Initial conditions
initState <- c(S = N - 1, I = 1, R = 0)

# Choose four values of v
#v.values <- c(1.0896, 1.37034, 1.7735, 2.7073)#MEAN
v.values <- c(0.98452, 1.11517, 1.31708, 1.82484)#specifc state

# Run simulations for each v
results <- list()

for (v in v.values) {
  params <- c(beta = beta, gamma = gamma, p = 1 + v^2, N = N)
  out <- ode(initState, t, SIRode, params)
  out <- as.data.frame(out)
  out$v <- v
  results[[as.character(v)]] <- out
}

# Combine results
all_out <- do.call(rbind, results)

# Plot
library(ggplot2)

ggplot(all_out, aes(x = time, y = I, color = factor(v))) +
  geom_line(size = 1) +
  #labs(title = "Deterministic SIR Epidemic Model",
  labs(title = "Deterministic SIR Epidemic Model Oyo",#specific states
       x = "Time",
       y = "Observed Infected",
       color = "v value") +
  theme_minimal(base_size = 14)#model1

##v and p_reportejoin
# Load packages
library(deSolve)
library(ggplot2)

# Parameters
N <- 4000000         # total population
gamma <- 0.2         # recovery rate
beta <- 0.41         # transmission rate
t <- seq(0, 200, by = 1)   # time steps

# Define the SIR model
SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS.dt <- -beta * (S/N)^p * I
    dI.dt <-  beta * (S/N)^p * I - gamma * I
    dR.dt <-  gamma * I
    return(list(c(dS.dt, dI.dt, dR.dt)))
  })
}

# Initial conditions
initState <- c(S = N - 1, I = 1, R = 0)

# Choose values for v and p_reported
v.values <- c(1.0896, 1.37034, 1.7735, 2.7073)
p.reported.values <- c(0.001, 0.005, 0.01, 0.02)

# Run simulations
results <- list()

for (v in v.values) {
  for (p_reported in p.reported.values) {
    params <- c(beta = beta, gamma = gamma, p = 1 + v^2, N = N)
    out <- ode(initState, t, SIRode, params)
    out <- as.data.frame(out)
    out$v <- v
    out$p_reported <- p_reported
    out$Observed <- out$I * p_reported   # scale infections by reporting fraction
    results[[paste(v, p_reported, sep = "_")]] <- out
  }
}

# Combine all results
all_out <- do.call(rbind, results)

# Plot observed cases
ggplot(all_out, aes(x = time, y = Observed, 
                    color = factor(p_reported), 
                    linetype = factor(v))) +
  geom_line(size = 1) +
  labs(title = "Observed Infected Cases (SIR Model)",
       x = "Time",
       y = "Observed Cases",
       color = "p_reported",
       linetype = "v value") +
  theme_minimal(base_size = 14)
#end 

##Model 2
# Load packages
library(deSolve)
library(ggplot2)

# Parameters
N <- 4000000         # total population
gamma <- 0.2         # recovery rate
R0 <- 2.5            # pick a fixed R0 (you can change later)
t <- seq(0, 200, by = 1)   # time steps

# Define the SIR model with R0 and pN only
SIRode <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma   # recompute beta from R0 and gamma
    dS <- -beta * I * (S / (pN * N))
    dI <- beta * I * (S / (pN * N)) - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
initState <- c(S = N - 1, I = 1, R = 0)

# Values of pN
pN.values <- c(0.72824149, 0.72739871, 0.72750674, 0.72750)

# Run simulations for each pN
results <- list()

for (pN in pN.values) {
  params <- c(R0 = R0, gamma = gamma, pN = pN, N = N)
  out <- ode(initState, t, SIRode, params)
  out <- as.data.frame(out)
  out$pN <- pN
  results[[as.character(pN)]] <- out
}

# Combine results
all_out <- do.call(rbind, results)

# Plot: Infected curves for each pN
ggplot(all_out, aes(x = time, y = I, color = factor(pN))) +
  geom_line(size = 1) +
  labs(title = "Deterministic SIR Epidemic Model (different pN values)",
       x = "Time",
       y = "Infected",
       color = "pN value") +
  theme_minimal(base_size = 14)# end model 2
##

##4 R0 delince
library(deSolve)
library(ggplot2)

# Parameters
N <- 4000000
gamma <- 0.2
beta <- 0.41
R0 <- beta / gamma   # baseline R0
t <- seq(0, 200, by = 1)

# SIR model with v
SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * (S/N)^p * I
    dI <-  beta * (S/N)^p * I - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
initState <- c(S = N - 1, I = 1, R = 0)

# v values
v.values <- c(1.0896, 1.37034, 1.7735, 2.7073)

results <- list()

for (v in v.values) {
  params <- c(beta = beta, gamma = gamma, p = 1 + v^2, N = N)
  out <- ode(initState, t, SIRode, params)
  out <- as.data.frame(out)
  
  # Compute Re(t) = R0 * (S/N)^(1+v^2)
  out$Re <- R0 * (out$S / N)^(1 + v^2)
  out$v <- v
  
  results[[as.character(v)]] <- out
}

all_out <- do.call(rbind, results)

# Plot
ggplot(all_out, aes(x = time, y = Re, color = factor(v))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Effective Reproduction Number over Time",
       x = "Time",
       y = expression(R[e](t)),
       color = "v value") +
  theme_minimal(base_size = 14)#end



##Specific states~EDO R0 delince
library(deSolve)
library(ggplot2)

# Parameters
N <- 3192893
gamma <- 0.2
beta <- 0.41
R0 <- beta / gamma   # baseline R0
t <- seq(0, 200, by = 1)

# SIR model with v
SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * (S/N)^p * I
    dI <-  beta * (S/N)^p * I - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
initState <- c(S = N - 1, I = 1, R = 0)

# v values
v.values <- c(0.96801, 1.26053, 1.59433, 2.35342)

results <- list()

for (v in v.values) {
  params <- c(beta = beta, gamma = gamma, p = 1 + v^2, N = N)
  out <- ode(initState, t, SIRode, params)
  out <- as.data.frame(out)
  
  # Compute Re(t) = R0 * (S/N)^(1+v^2)
  out$Re <- R0 * (out$S / N)^(1 + v^2)
  out$v <- v
  
  results[[as.character(v)]] <- out
}

all_out <- do.call(rbind, results)

# Plot
ggplot(all_out, aes(x = time, y = Re, color = factor(v))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Effective Reproduction Number over Time for Kwara",
       x = "Time",
       y = expression(R[e](t)),
       color = "v value") +
  theme_minimal(base_size = 14)#end


###Ro decline for model 2
## Decline in Re(t) with pN formulation
library(deSolve)
library(ggplot2)

# Parameters
N <- 3192893
gamma <- 0.2
R0 <- 2.05   # choose baseline R0 (you can adjust)
t <- seq(0, 200, by = 1)

# Define the SIR model with R0 and pN
SIRode <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma   # recompute beta from R0 and gamma
    dS <- -beta * I * (S / (pN * N))
    dI <-  beta * I * (S / (pN * N)) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
initState <- c(S = N - 1, I = 1, R = 0)

# Values of pN
pN.values <- c(0.73482, 0.73205, 0.73249, 0.73347)

results <- list()

for (pN in pN.values) {
  params <- c(gamma = gamma, R0 = R0, pN = pN, N = N)
  out <- ode(initState, t, SIRode, params)
  out <- as.data.frame(out)
  
  # Compute Re(t) = R0 * (S / (pN * N))
  out$Re <- R0 * (out$S / (pN * N))
  out$pN <- pN
  
  results[[as.character(pN)]] <- out
}

all_out <- do.call(rbind, results)

# Plot
ggplot(all_out, aes(x = time, y = Re, color = factor(pN))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Delcine in Effective Reproduction Number for Kwara",
       x = "Time",
       y = expression(R[e](t)),
       color = "pN value") +
  theme_minimal(base_size = 14)#end



#####INCLUDING OBSERVED CASES TO THE PLOT FOR MODEL 1
library(deSolve)
library(ggplot2)

# Parameters
N <- 4000000         
gamma <- 0.2 
beta <- 0.41
t <- seq(0, 200, by = 1)   

SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS.dt <- -beta * (S / N)^p * I
    dI.dt <-  beta * (S / N)^p * I - gamma * I
    dR.dt <-  gamma * I
    incidence <- -dS.dt               # daily new cases
    cum_cases <- I + R                # cumulative infected
    return(list(c(dS.dt, dI.dt, dR.dt),
                incidence = incidence,
                cum_cases = cum_cases))
  })
}

initState <- c(S = N - 1, I = 1, R = 0)

v.values <- c(0.96801, 1.26053, 1.59433, 2.35342)  # specific states and scenarios Oyo
#simu
results <- list()

for (v in v.values) {
  params <- c(beta = beta, gamma = gamma, p = 1 + v^2, N = N)
  out <- ode(initState, t, SIRode, params)
  out <- as.data.frame(out)
  out$v <- v
  results[[as.character(v)]] <- out
}

# results into a data frame
all_out <- do.call(rbind, results)

# Plot 1:Infectious (I) 
ggplot(all_out, aes(x = time, y = I, color = factor(v))) +
  geom_line(size = 1) +
  labs(title = "Infectious (I) - Deterministic model Kwara",
       x = "Time",
       y = "Number of Infectious Individuals",
       color = "v value") +
  theme_minimal(base_size = 14)

# Plot 2: Daily New Cases (Incidence) 
ggplot(all_out, aes(x = time, y = incidence, color = factor(v))) +
  geom_line(size = 1) +
  labs(title = "Daily New Cases (Incidence) - Det model Kwara",
       x = "Time",
       y = "New Cases per Day",
       color = "v value") +
  theme_minimal(base_size = 14)

# Plot 3: Cumulative Number of Cases 
ggplot(all_out, aes(x = time, y = cum_cases, color = factor(v))) +
  geom_line(size = 1) +
  labs(title = "Cumu Number of Cases - Det model kwara",
       x = "Time",
       y = "Cumulative Infections",
       color = "v value") +
  theme_minimal(base_size = 14)#end model 1 including incidence  (dialy cases)

