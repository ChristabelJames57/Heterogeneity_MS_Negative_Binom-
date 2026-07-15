#library(posterior distribution)
library(bayesplot)

#save in
# Extract all posterior draws
#draws <- stan_out$draws()

# Save them for reproducibility
#saveRDS(draws, "./Output/draws_model11_heterogeneity_NB.rds")
#saveRDS(selected_states, "./Output/selected_states.rds")

#read back in
draws <- readRDS("./Output/draws_model11_heterogeneity_NB.rds")
selected_states <- readRDS("./Output/selected_states.rds")


draws <- stan_out$draws()
#plot of parameters
mcmc_dens(
  draws,
  pars = c("MeanV",
           "meanbeta",
           "meanp_reported",
           "meanR0",
           "meanI0")
)

#for Initial conditions
mcmc_dens(
  draws,
  pars = paste0("I0[",1:n_states,"]")
)

#option2 Histogram of postrior dsn
mcmc_hist(
  draws,
  pars = c("MeanV",
           "meanbeta",
           "meanR0",
           "meanp_reported",
           "meanI0")
)

#or
mcmc_areas(
  draws,
  pars = c("MeanV",
           "meanbeta",
           "meanR0",
           "meanp_reported",
           "meanI0"),
  prob = 0.95
)



##another
library(tidybayes)
library(ggplot2)

posterior <- stan_out %>%
  spread_draws(MeanV,
               meanbeta,
               meanR0,
               meanp_reported,
               meanI0)

ggplot(posterior, aes(x = MeanV)) +
  geom_density(fill="skyblue", alpha=.5) +
  theme_bw()

ggplot(posterior, aes(x = meanI0)) +
  geom_density(fill="skyblue", alpha=.5) +
  theme_bw()


##save posterior draws
posterior_df <- as_draws_df(draws)

head(posterior_df)
write.csv(
  posterior_df,
  "PosteriorSamples.csv",
  row.names = FALSE
)

##now new 
############################################################
## POSTERIOR DISTRIBUTIONS FOR ALL PARAMETERS
############################################################

library(cmdstanr)
library(bayesplot)
library(posterior)
library(ggplot2)

#-----------------------------------------------------------
# Extract posterior draws
#-----------------------------------------------------------

draws <- stan_out$draws()

# Convert to posterior dataframe
posterior_df <- as_draws_df(draws)

# Save posterior samples
write.csv(posterior_df,
          "PosteriorSamples.csv",
          row.names = FALSE)

#-----------------------------------------------------------
# List all parameters
#-----------------------------------------------------------

all_pars <- variables(draws)

print(all_pars)

############################################################
## 1. Density plots for ALL parameters
############################################################

pdf("Posterior_Densities.pdf", width=8, height=6)

for(par in all_pars){
  
  print(
    mcmc_dens(
      draws,
      pars = par
    ) +
      ggtitle(paste("Posterior Distribution of", par))
  )
  
}

dev.off()

############################################################
## 2. Histograms
############################################################

pdf("Posterior_Histograms.pdf", width=8, height=6)

for(par in all_pars){
  
  print(
    mcmc_hist(
      draws,
      pars = par
    ) +
      ggtitle(paste("Posterior Histogram of", par))
  )
  
}

dev.off()

############################################################
## 3. Area plots (95% Credible Intervals)
############################################################

pdf("Posterior_Areas.pdf", width=8, height=6)

for(par in all_pars){
  
  print(
    mcmc_areas(
      draws,
      pars = par,
      prob = 0.95
    ) +
      ggtitle(paste("Posterior Distribution of", par))
  )
  
}

dev.off()

############################################################
## 4. Trace plots (Convergence)
############################################################

pdf("Trace_Plots.pdf", width=8, height=6)

for(par in all_pars){
  
  print(
    mcmc_trace(
      draws,
      pars = par
    ) +
      ggtitle(paste("Trace Plot:", par))
  )
  
}

dev.off()

############################################################
## 5. Numerical posterior summary
############################################################

posterior_summary <- summarise_draws(draws)

print(posterior_summary)

write.csv(
  posterior_summary,
  "Posterior_Summary.csv",
  row.names = FALSE
)

############################################################
## 6. Posterior distributions of key parameters together
############################################################

key_pars <- c(
  "MeanV",
  "meanbeta",
  "meanR0",
  "meanp_reported",
  "meanI0",
  "theta_ss",
  "slab_sd",
  "R0"
)

key_pars <- key_pars[key_pars %in% all_pars]

mcmc_areas(
  draws,
  pars = key_pars,
  prob = 0.95
)

############################################################
## 7. State-specific beta
############################################################

beta_pars <- grep("^beta\\[", all_pars, value=TRUE)

if(length(beta_pars)>0){
  
  mcmc_areas(
    draws,
    pars = beta_pars,
    prob = 0.95
  )
  
}

############################################################
## 8. State-specific v
############################################################

v_pars <- grep("^v\\[", all_pars, value=TRUE)

if(length(v_pars)>0){
  
  mcmc_areas(
    draws,
    pars = v_pars,
    prob = 0.95
  )
  
}

############################################################
## 9. State-specific I0
############################################################

I0_pars <- grep("^I0\\[", all_pars, value=TRUE)

if(length(I0_pars)>0){
  
  mcmc_areas(
    draws,
    pars = I0_pars,
    prob = 0.95
  )
  
}

############################################################
## 10. Posterior means and credible intervals
############################################################

posterior_summary %>%
  select(variable,
         mean,
         median,
         sd,
         q5,
         q95,
         rhat,
         ess_bulk,
         ess_tail)

############################################################
## Finished
############################################################