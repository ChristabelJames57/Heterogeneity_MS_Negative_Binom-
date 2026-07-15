library(tidyverse)
library(posterior)
library(grid)
library(deSolve)
library(scales)
library(patchwork)
rm(list = ls()) # Clear workspace

##saved
# Save posterior draws for model 1
#draws_df_m11 <- as_draws_df(stan_out$draws("pred_cases"))
#saveRDS(draws_df_m11, "./Output/draws_model11_heterogeneity_NB.rds")#mhet
#draws_df_m33 <- as_draws_df(stan_out$draws("pred_cases"))
#saveRDS(draws_df_m33, "./Output/draws_model33_heterogeneity_NB.rds")#mhom
#saveRDS(draws_df_m44, "./Output/draws_model44_heterogeneity_NB.rds")#mhom-het
#saveRDS(selected_states, "./Output/selected_states.rds")


# Figure_6 reported versus predicted for all three models #####

draws_df_m11 <- readRDS("./Output/draws_model11_heterogeneity_NB.rds")
draws_df_m33 <- readRDS("./Output/draws_model33_noheterogeneity_NB.rds")
draws_df_m44 <- readRDS("./Output/draws_model44_PNheterogeneity_NB.rds")

cases_matrix    <- readRDS("./Output/cases_matrix.rds")
selected_states <- readRDS("./Output/selected_states.rds")

#   t_last 
end_days <- c(
  Ondo = 152, Edo = 189, Osun = 127, Oyo = 174,
  Kano = 144, Kwara = 154, Delta = 140,
  Ebonyi = 119, Enugu = 128, Kaduna = 200, Bayelsa = 89,
  Gombe = 143, Rivers = 157,
  Bauchi = 129, `Cross River` = 89
)

t_last <- end_days[selected_states]  # correct order

# SUMMARISE FUNCTION

summarise_predictions <- function(draws_df, t_last, model_name) {
  
  n_states <- length(t_last)
  
  t_last_df <- tibble(
    state = 1:n_states,
    t_last = t_last
  )
  
  draws_df %>%
    pivot_longer(
      cols = starts_with("pred_cases["),
      names_to = "var",
      values_to = "value"
    ) %>%
    mutate(
      var = str_remove_all(var, "pred_cases\\[|\\]"),
      state = as.integer(str_extract(var, "^[0-9]+")),
      day   = as.integer(str_extract(var, "(?<=,)[0-9]+"))
    ) %>%
    left_join(t_last_df, by = "state") %>%
    filter(day <= t_last) %>%  
    group_by(state, day) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(model = model_name)
}

# SUMMARISE ALL MODELS

summary_m11 <- summarise_predictions(draws_df_m11, t_last, "Heterogeneity")
summary_m33 <- summarise_predictions(draws_df_m33, t_last, "No heterogeneity")
summary_m44 <- summarise_predictions(draws_df_m44, t_last, "Mhom-het")

summary_all <- bind_rows(summary_m11, summary_m33, summary_m44)

# OBSERVED DATA (FIXED)

cases_subset <- as.matrix(cases_matrix[1:max(t_last), ])

obs_df <- purrr::map_dfr(1:length(t_last), function(s) {
  tibble(
    state = s,
    day = 1:t_last[s],  
    Observed = cases_subset[1:t_last[s], s]
  )
})

# MERGE

plot_df <- summary_all %>%
  left_join(obs_df, by = c("state", "day"))

# Labels
state_labels <- setNames(
  selected_states,
  as.character(1:length(selected_states))
)

# PLOT

ggplot(plot_df, aes(x = day)) +
  geom_point(
    aes(y = Observed),
    color = "black",
    size = 0.8
  ) +
  geom_line(
    aes(y = median, color = model, linetype = model),
    linewidth = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Heterogeneity"     = "red",
      "No heterogeneity"  = "blue",
      "Mhom-het"          = "darkgreen"
    )
  ) +
  facet_wrap(
    ~state,
    labeller = as_labeller(state_labels),
    scales = "free" 
  ) +
  
  labs(
    x = "Day",
    y = "Reported cases",
    color = "Model"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save
ggsave("Output/Figure/Figure_6.pdf",
       width = 12, height = 10,
       units = "in", bg = "white")


#Figure_7 Deterministic outbreak comparison using the mean states #####
#mean states and four models


### MODEL: SIR with heterogeneity
SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    incidence <- beta * (S / N)^p * I^q
    
    dS.dt <- -incidence
    dI.dt <-  incidence - gamma * I
    dR.dt <-  gamma * I
    
    return(list(
      c(dS.dt, dI.dt, dR.dt),
      reported = incidence
    ))
  })
}

### FIXED PARAMETERS
gamma <- 0.2
q <- 1
times <- seq(1, 330, by = 1)

### POPULATION
N  <- 5315779
PN <- 0.508

##### 1. Mhom (p = 1)

beta1 <- 0.272
I01   <- 1052.8
S01   <- N - I01

params0 <- list(
  beta = beta1,
  gamma = gamma,
  p = 1,
  q = q,
  N = N
)

out0 <- as.data.frame(
  ode(c(S = S01, I = I01, R = 0),
      times, SIRode, params0)
)

out0$scenario <- "Mhom (p = 1)"


##### 2. Mhet

v1 <- 2.561

params1 <- list(
  beta = beta1,
  gamma = gamma,
  p = 1 + v1^2,
  q = q,
  N = N
)

out1 <- as.data.frame(
  ode(c(S = S01, I = I01, R = 0),
      times, SIRode, params1)
)

out1$scenario <- "Mhet"


##### 3. Mhom-het

beta2 <- 0.272
v2    <- 2.001
I02   <- 728.3

N_eff <- N * PN
S02   <- N_eff - I02

params2 <- list(
  beta = beta2,
  gamma = gamma,
  p = 1 + v2^2,
  q = q,
  N = N_eff
)

out2 <- as.data.frame(
  ode(c(S = S02, I = I02, R = 0),
      times, SIRode, params2)
)

out2$scenario <- "Mhom-het"


##### 4. Mhom (p = 0.508)

S03 <- N_eff - I01

params3 <- list(
  beta = beta1,
  gamma = gamma,
  p = 1,
  q = q,
  N = N_eff
)

out3 <- as.data.frame(
  ode(c(S = S03, I = I01, R = 0),
      times, SIRode, params3)
)

out3$scenario <- "Mhom (p = 0.508)"


##### COMBINE ALL

df_plot <- bind_rows(out0, out1, out2, out3)

#nforce legend order
df_plot$scenario <- factor(
  df_plot$scenario,
  levels = c(
    "Mhom-het",
    "Mhet",
    "Mhom (p = 0.508)",
    "Mhom (p = 1)"
  )
)

##### COLOURS

scenario_colours <- c(
  "Mhom-het"         = "darkgreen",
  "Mhet"             = "red",
  "Mhom (p = 0.508)" = "purple",
  "Mhom (p = 1)"     = "blue"
)

###### PLOT

ggplot(df_plot,
       aes(x = time, y = reported,
           colour = scenario)) +
  
  geom_line(linewidth = 1.2) +
  
  scale_colour_manual(values = scenario_colours) +
  
  labs(
    x = "Time",
    y = "Reported Infected",
    colour = "Scenario"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  )
# Save
ggsave("Output/Figure/Figure_7.pdf", width = 7, height = 5, bg = "white")
# END Figure_7


#Figure_8 Gamma dsn with a mean of 1 and cv=2.561 and Herd Immunity threshold##### 
gamma <- 0.2

beta_vec <- c(
  0.3009, 0.2935, 0.2676, 0.2896, 0.2367,
  0.2656, 0.2914, 0.2689, 0.2565, 0.2506
)

v_vec <- c(
  3.0728, 2.5497, 1.7406, 3.1850, 1.8934,
  2.2834, 2.8745, 3.3145, 1.7874, 2.8323
)

R0_vec <- beta_vec / gamma

state_df <- data.frame(
  state = factor(1:10),
  R0    = R0_vec,
  v_hat = v_vec
)

v_seq <- seq(0, 4, length.out = 400)

hit_df <- expand.grid(
  v     = v_seq,
  state = state_df$state
) %>%
  left_join(state_df, by = "state") %>%
  mutate(
    HIT = 1 - R0^(-1 / (1 + v^2))
  )

hit_points <- state_df %>%
  mutate(
    HIT = 1 - R0^(-1 / (1 + v_hat^2))
  )

CV <- 2.561
shape <- 1 / (CV^2)
scale <- 1 / shape

x <- seq(0, 20, length.out = 1000)
df_gamma <- data.frame(
  x = x,
  y = dgamma(x, shape = shape, scale = scale)
)

p_inset <- ggplot(df_gamma, aes(x = x, y = y)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Susceptibility",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Convert inset to grob
g_inset <- ggplotGrob(p_inset)

# == FINAL PLOT WITH INSET
p_final <- ggplot(hit_df, aes(
  x = v,
  y = HIT,
  group = state
)) +
  
  geom_line(colour = "black", linewidth = 0.9) +
  
  geom_point(
    data = hit_points,
    aes(x = v_hat, y = HIT),
    size = 3,
    colour = "red"
  ) +
  
  labs(
    x = expression("Coefficient of variation (" * v * ")"),
    y = "Herd immunity threshold"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  # Inset
  annotation_custom(
    grob = g_inset,
    xmin = 2.5, xmax = 3.9,
    ymin= 0.09, ymax= 0.33
  )

p_final

# Save
ggsave("Output/Figure/Figure_8.pdf", width = 7, height = 5, bg = "white")
#end figure_8


#Figure_9 DECLINE IN EFFECTIVE REPRODUCTION NUMBER (mean states) four models (panels)#####
library(scales)
library(patchwork)

###### SIR MODEL ######

SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    incidence <- beta * (S / N)^p * I^q
    
    dS <- -incidence
    dI <-  incidence - gamma * I
    dR <-  gamma * I
    
    list(c(dS, dI, dR))
  })
}

#### PARAMETERS ####

gamma <- 0.2
q <- 1
times <- seq(1, 330, by = 1)

###### MEAN POPULATION ######

N  <- 5315779
PN <- 0.508

##### 1. Mhom (p = 1) #####

beta  <- 0.272  
v     <- 2.561 
I0    <- 1052.8 

S0 <- N - I0
init0 <- c(S = S0, I = I0, R = 0)

R0_basic0 <- beta / gamma
p0 <- 1

out0 <- as.data.frame(
  ode(init0, times, SIRode,
      parms = list(beta = beta,
                   gamma = gamma,
                   p = p0,
                   q = q,
                   N = N))
)

out0$model <- "Mhom (p = 1)"
out0$Rt <- R0_basic0 * (out0$S / N)^p0


##### 2. Mhet #####

p1 <- 1 + v^2

out1 <- as.data.frame(
  ode(init0, times, SIRode,
      parms = list(beta = beta,
                   gamma = gamma,
                   p = p1,
                   q = q,
                   N = N))
)

out1$model <- "Mhet"
out1$Rt <- R0_basic0 * (out1$S / N)^p1


##### 3. Mhom-het #####

beta2 <- 0.272 
v2    <- 2.001 
I02   <- 728.3 

N_eff2 <- N * PN
S02    <- N_eff2 - I02

init2 <- c(S = S02, I = I02, R = 0)

R0_basic2 <- beta2 / gamma
p2 <- 1 + v2^2

out2 <- as.data.frame(
  ode(init2, times, SIRode,
      parms = list(beta = beta2,
                   gamma = gamma,
                   p = p2,
                   q = q,
                   N = N_eff2))
)

out2$model <- "Mhom-het"
out2$Rt <- R0_basic2 * (out2$S / N_eff2)^p2


##### 4. Mhom (p = 0.508) #####

N_eff0 <- N * PN
S03    <- N_eff0 - I0

init3 <- c(S = S03, I = I0, R = 0)

out3 <- as.data.frame(
  ode(init3, times, SIRode,
      parms = list(beta = beta,
                   gamma = gamma,
                   p = 1,
                   q = q,
                   N = N_eff0))
)

out3$model <- "Mhom (p = 0.508)"
out3$Rt <- R0_basic0 * (out3$S / N_eff0)


##### COMBINE #####

df_all <- bind_rows(out0, out1, out2, out3)

df_all$model <- factor(
  df_all$model,
  levels = c(
    "Mhom-het",
    "Mhet",
    "Mhom (p = 0.508)",
    "Mhom (p = 1)"
  )
)

##### COMPUTE TIME TO Rt = 1 #####

time_to_R1 <- df_all %>%
  arrange(model, time) %>% ##like sorting
  group_by(model) %>% ##each model stands seperately
  summarise(
    t_R1 = {     #creating one number per model
      idx <- which(Rt <= 1)[1]   ##time to which epidemic is declining first outomce is 1
      
      if (is.na(idx) || idx == 1) { #Rt never drops below 1 and already  less than  or equal to 1 at start
        NA   ##both returns NA
      } else {
        t1 <- time[idx - 1]  ##Gets the two points arounnd
        t2 <- time[idx]
        r1 <- Rt[idx - 1]
        r2 <- Rt[idx]
        
        t1 + (1 - r1) * (t2 - t1) / (r2 - r1)
      }
    }
  )

print(time_to_R1)

##### COLOURS #####

scenario_colours <- c(
  "Mhom-het"         = "darkgreen",
  "Mhet"             = "red",
  "Mhom (p = 0.508)" = "purple",
  "Mhom (p = 1)"     = "blue"
)

######## PANEL A ########

p1 <- ggplot(df_all, aes(x = R, y = Rt, colour = model)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_continuous(labels = comma) +
  scale_colour_manual(values = scenario_colours) +
  labs(
    x = "Recovered",
    y = expression(R[t]),
    colour = "Scenario",
    title = "A"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

######## PANEL B (UPDATED) ########

p2 <- ggplot(df_all, aes(x = time, y = Rt, colour = model)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  # vertical lines at Rt = 1 crossing
  geom_vline(
    data = time_to_R1,
    aes(xintercept = t_R1, colour = model),
    linetype = "dotted",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(values = scenario_colours) +
  labs(
    x = "Time (days)",
    y = expression(R[t]),
    colour = "Scenario",
    title = "B"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )
p1 / p2
# Save
ggsave("Output/Figure/Figure_9.pdf", width = 7, height = 5, bg = "white")
#end Figure_9


