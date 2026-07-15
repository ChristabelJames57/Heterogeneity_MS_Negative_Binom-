#Mhet model using Negative Binom observation error##### 
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

#Mhet model######

set.seed(1275) #real
#set.seed(1275) simu
total_number_of_time_series=19#for simulation 15NL and 3L must bge highr than noof outbvreak
MAX_TIME<-189 ## To test real data
N=matrix(data=0, nrow = total_number_of_time_series)
tru_incidence_matrix <- obs_incidence_matrix<- matrix(data=0, nrow = MAX_TIME , ncol = total_number_of_time_series)#all in a matrix with that row n column
Incidence=matrix(data=0, nrow = MAX_TIME,ncol=total_number_of_time_series)#so we have day row n outbreak col
obs_incidence_matrix2 <- matrix(0, nrow = MAX_TIME, ncol = total_number_of_time_series)


SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma  # recompute beta from R0 and gamma
    dS <- -beta * I * (S / N[i])^(1 + v^2)
    dI <-  beta * I * (S / N[i])^(1 + v^2) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}


cmdstanr_max_rows=15
max_rows = 15
#setwd("~/Desktop/Dans Folder 2023/People/Christabel James/Stan code/Attempt3")

The_mode="S"#change only here from  S to R N chge prior for rp
if (The_mode=="R"){
  df_N <- read.csv("./Data/Nigeria17.csv") # Columns: Day, Edo,Lagos, Ondo, Ogun, Osun, Oyo, 
  populations <- c(Ondo = 4671695, Edo = 4235595, Osun = 4705589, Ogun = 5217716, Oyo = 7840864, Lagos = 12550598, Kano = 13076892, Kwara = 3192893, Delta = 5663362, Ebonyi = 2880383, Enugu = 4411119, Kaduna = 8252366, Bayelsa = 2277961, Benue = 5741815, Fct = 3564126, Gombe = 3256962, Rivers = 7303924, Bauchi = 6537314, `Cross River` = 3866269)
  end_days<-c(Ondo = 152, Edo = 189, Osun = 127, Ogun = 197, Oyo = 174, Lagos = 228,Kano =144, Kwara = 154, Delta = 140, Ebonyi = 119, Enugu = 128, Kaduna = 200, Bayelsa =89, Benue =86, Fct= 225, Gombe = 143, Rivers= 157, Bauchi = 129, `Cross River` =89)#real data 15 states

  
  # Select States (ten)
  selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara","Ondo", "Osun", "Oyo", "Rivers")#arrge alphabe 13 no Benue
  
  # Prepare Data 
  n_states <- length(selected_states)
  #n_days <- length(df_N$Day)
  t_last <- end_days[selected_states]#end days for selcted states alone
  
  cases_matrix <- df_N %>%#matrix of case data
    select(all_of(selected_states)) %>%#selects the column for specified state
    as.matrix(nrow=t_last,ncol=n_states)#converts to matrix
  
  
  N_vector <- populations[selected_states]
  #end real data
  
} else{#begin simulated data

  beta <- 0.4#
  gamma <- 0.2
  #pN=0.4
  v <-2 #Fixing homogenous
  R0 <- beta / gamma
  reporting_prob <-  0.015#0.015#mean pr testing 1/pr-
  I0 <- rpois(1,1/ reporting_prob)#DIFFERENT Initial conditions
 
  sim_r <- 2.1
  noof_outbreaks <- 10 # using 10 simulated outbreaks as u
  n_states <- noof_outbreaks
  selected_states <- paste0("SimState_", 1:noof_outbreaks)
  
  # Delta, Ebonyi, Edo, Enugu, Kaduna, Kwara, Ondo, Osun, Oyo, Rivers
  pop_order <- c(
    5663362, # Delta
    2880383, # Ebonyi
    4235595, # Edo
    4411119, # Enugu
    8252366, # Kaduna
    3192893, # Kwara
    4671695, # Ondo
    4705589, # Osun
    7840864, # Oyo
    7303924  # Rivers
  )
  
  # pre-allocate vectors / matrices
  N_vector <- integer(noof_outbreaks)        # hold integer populations passed to Stan
  t_last <- integer(noof_outbreaks)

  if (!is.numeric(N)) N <- as.numeric(N)
  
  if (length(N) < noof_outbreaks) N <- numeric(noof_outbreaks)  # ensure N has at least noof_outbreaks values
  
  for (i in 1:noof_outbreaks) {
    N[i] <- pop_order[i]
    S0 <- N[i] - I0
    R0_init <- 0
    y_init <- c(S = S0, I = I0, R = R0_init)
    
    times <- seq(0, MAX_TIME, by = 1)
    
    parameters <- c(beta = beta, gamma = gamma, v = v)
    
    out <- ode(y = y_init, times = times, func = SIRode, parms = parameters)
    out_df <- as.data.frame(out)
    
    incidence <- -diff(out_df$S) # new infections
    incidence[incidence < 0] <- 0
    
    # True process
    tru_incidence_matrix[, i] <- incidence
    
    # Observed process (Neg-Bin)
    obs_incidence_matrix[, i] <- rnbinom(MAX_TIME, mu = reporting_prob * incidence + 1e-6, size = sim_r)# nEGATIVEBINOM
    #obs_incidence_matrix[, i] <- rbinom(MAX_TIME, size = (trunc(incidence)+1), prob = reporting_prob)#binom epid used
    obs_incidence_matrix2[, i] <- obs_incidence_matrix[, i]
    
    N_vector[i] <- as.integer(N[i])    # store population and t_last
    t_last[i] <- MAX_TIME
  }
  
  cases_matrix <- obs_incidence_matrix[, 1:noof_outbreaks] # matrix of cases
  

  
}#end simulated


#N = as.integer(N_vector)# ensure integers#SImulation added converts N vectors to integers after tlast

data_sir <- list(max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,N = N_vector,cases = cases_matrix[1:max(t_last), ], spike_sd = 1e-3)#pas as input stan no phi and REAL estim IO

#NCovid <- cmdstan_model("NCovid_Multi_Normal_FixvPr_Real_Estim Pr_RemoveIO_Global1.stan")#state speicif _GOODFEBRUARY2026#NCovid <- cmdstan_model("NCovid_Multi_nb_rpt_v_w.stan")# compile model suscep no random effc
NCovid <- cmdstan_model("./Scripts/Mhet_model_NB.stan")#-ve binom senstivity

iters=2000 #4000 for simu
stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4, seed = 0)
#stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4, seed = 0, adapt_delta = 0.99, max_treedepth = 15)# non convergen

print(stan_out, variables = c("MeanV", "meanbeta", "meanp_reported","meanR0", "meanI0", "theta_ss", "I0", "v", "beta", "phi"), digits = 8,max_rows = 70)


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
state_labels <- setNames(selected_states, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))#replac 1234 with abc for four states

# Plot using facets
ggplot(plot_df, aes(x = day)) +
  geom_point(aes(y = Observed), color = "black", size = 0.8) +
  geom_line(aes(y = median), color = "blue", size = 1) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.15)+#to lighten the colour for bette visbility
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +##CI
  #facet_wrap(~state, scales = "free_y") +
  facet_wrap(~state, labeller = as_labeller(state_labels), scales = "free_y") +#FORlabeling states
  labs(x = "Day", y = "Cases", title = "Observed vs Predicted Cases by State") +
  theme_minimal(base_size = 12)#end







###To save
library(posterior)

draws_df_m11 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m11, file = "draws_model11_heterogeneity_NB.rds")

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