#Figure 1: Incidence data for the first wave #####
library(tidyverse)
library(plotly)
library(ggforce)
library(zoo)
library(sf)
library(posterior)
library(grid)
library(deSolve)
library(scales)
library(patchwork)

rm(list = ls()) # Clear workspace

##Figure_1 
#Data_wave1 <- read_csv("./Data/wave1_data.csv")
Data_wave12 <- read_csv("./Data/wave12_data.csv")#new

###Data has day, month and year can create a proper date column
Data_wave12 <- Data_wave12 %>%
  mutate(Date = make_date(year, month, day)) 
names(Data_wave12)


# Aggregate data by date

D_State <- Data_wave12 %>%
  group_by(Date,`State of Residence`) %>% #can do multiple group by of variables
  summarize(Case_Count = n()) %>% 
  na.omit(`State of Residence`) #To remove the NA in state

#To check if any NA
any(is.na(D_State))# in the entire data frame of D_state. if NA will return true if no will return False
any(is.na(D_State$date))# to check for NA in date column
any(is.na(D_State$Case_Count))#check for NA in the case_count column


#Aggregate Cases by date and state
cases_by_state <- Data_wave12 %>% ##Shows the data frame to use for plotting
  group_by(State = `State of Residence`, Date) %>%
  summarise(Case_Count = n(), .groups = "drop")

###clean
cases_by_state_clean <- na.omit(cases_by_state)#remove NA before ploting


# Define a vector with the specific states to plot
#selected_states <- c("Bayelsa", "Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara", "Ondo", "Oyo", "Rivers")  # Replace with the states you're interested in

# Define a vector with the specific states to plot Updated
selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara", "Ondo", "Oyo", "Osun", "Rivers")  # Replace with the states you're interested in


# Filter the data to include only the selected states
cases_by_state_filtered <- cases_by_state_clean %>%
  filter(State %in% selected_states)

plot <- ggplot(cases_by_state_filtered, aes(x = Date, y = Case_Count)) +
  geom_line(color = "blue", linewidth = 0.4) +
  
  # Facet layout similar to the image
  facet_wrap(~ State, scales = "free_y", ncol = 4) +
  
  labs(
    #title = "Time Series of Cases by State",
    x = "Date",
    y = "Number of Cases"
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    # Title styling
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    # Facet labels
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    # Axis text formatting
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 8),
    
    # Clean panel look
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

print(plot)
ggsave("Output/Figure/Figure_1_updated.pdf", width = 7, height = 5, bg = "white")
#END Figure 1

# Figure_2 MAP OF SELECTED STATES ######
# Load the shapefile using sf

nigeria_states <- st_read("Data/Shapefile/gadm41_NGA_1.shp")


#prepare data
# Check the column names
print(names(nigeria_states))

# Clean and format the state names
nigeria_states$name <- tools::toTitleCase(tolower(nigeria_states$NAME_1))


# Define a vector with the specific states to plot updated
selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara", "Ondo", "Oyo", "Osun", "Rivers")  # Replace with the states you're interested in


# Mark states as "Selected" or "Other"
nigeria_states <- nigeria_states %>%
  mutate(selected = ifelse(name %in% selected_states, "Selected", "Other"))

# Plot the map

ggplot(nigeria_states) +
  geom_sf(aes(fill = selected), color = "black", size = 0.2) +
  scale_fill_manual(values = c("Selected" = "skyblue", "Other" = "grey90")) +
  theme_minimal() +
  #labs(title = "Map of Nigeria Highlighting 17 Selected States",
  #fill = "State Category")
  labs(title = "",
       fill = "State Category")

ggsave("Output/Figure/Figure_2_Updated.pdf", width = 7, height = 5, bg = "white")
###END Fig.2 




#Figure_3 Recovery of R0 parameter#####

# Load data
df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_use.csv")

# Treat v as factor
df_R0 <- df_R0 %>%
  mutate(v = factor(v))

# Position dodge for consistent separation
pd <- position_dodge(width = 0.0008)

ggplot(df_R0, aes(x = p_r, y = MeanR0, colour = v, group = v)) +
  
  # Mean estimate
  geom_point(
    size = 3,
    position = pd
  ) +
  
  # 95% credible intervals
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0004,
    linewidth = 0.7,
    position = pd
  ) +
  
  # True R0 reference line
  geom_hline(
    yintercept = 2,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.8
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "blue",
      "1" = "green",
      "2" = "red",
      "3" = "black"
    )
  ) +
  
  labs(
    #x = expression(p[reported]),
    # y = expression(hat(R)[0])
    x = expression(italic(r)),
    y = expression(hat(R)[0])
  ) +
  
  scale_y_continuous(
    limits = c(1.80, 2.05),
    breaks = c(1.95, 2.00, 2.06)
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Keep clean look
    panel.grid = element_blank(), # remove all grid lines
    
    plot.margin = margin(10, 10, 10, 10)
  )

# Save
ggsave("Output/Figure/Figure_3.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_3.png", width = 7, height = 5, bg = "white", dpi=1000)
#end Figure_3



#using 0 to 2
#Figure_3 Recovery of R0 parameter#####

# Load data
#df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_use2.csv")#gam(1,10)for Normal and NB
#df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_Gamma.csv")#using gama (1,1)
df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_Gam01.csv")#using gama (0.1,0.1)

# Treat v as factor
df_R0 <- df_R0 %>%
  mutate(v = factor(v))

# Position dodge for consistent separation
pd <- position_dodge(width = 0.0008)

ggplot(df_R0, aes(x = p_r, y = MeanR0, colour = v, group = v)) +
  
  # Mean estimate
  geom_point(
    size = 3,
    position = pd
  ) +
  
  # 95% credible intervals
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0004,
    linewidth = 0.7,
    position = pd
  ) +
  
  # True R0 reference line
  geom_hline(
    yintercept = 2,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.8
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "blue",
      "1" = "green",
      "2" = "red"
    )
  ) +
  
  labs(
    #x = expression(p[reported]),
    # y = expression(hat(R)[0])
    x = expression(italic(r)),
    y = expression(hat(R)[0])
  ) +
  
  scale_y_continuous(
    #limits = c(1.80, 2.05),
    limits = c(1.80, 2.06),
    breaks = c(1.95, 2.00, 2.06)
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Keep clean look
    panel.grid = element_blank(), # remove all grid lines
    
    plot.margin = margin(10, 10, 10, 10)
  )

# Save
ggsave("Output/Figure/Figure_3.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_3.png", width = 7, height = 5, bg = "white", dpi=1000)
#end Figure_3



#Figure_4 Recovery of r parameter##### 

# Read data
#df <- read_csv("./Data/MeanPreported_summaryNEW_NB_use.csv")#0 to 3
df <- read_csv("./Data/MeanPreported_summaryNEW_NB_use2.csv")# 0 to 2
# Treat v as factor
df <- df %>%
  mutate(v = factor(v))

# Plot
ggplot(df, aes(x = p_r, y = Meanp_r, colour = v)) +
  
  # Mean estimates
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # 95% credible interval (horizontal dashes)
  geom_segment(
    aes(
      x = p_r - 0.00015, xend = p_r + 0.00015,
      y = LowerCI, yend = LowerCI
    ),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(
      x = p_r - 0.00015, xend = p_r + 0.00015,
      y = UpperCI, yend = UpperCI
    ),
    linewidth = 0.9
  ) +
  
  # True r reference lines
  geom_hline(
    aes(yintercept = p_r_true),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # Facet by ν (labels removed below)
  facet_wrap(~ v, ncol = 1) +
  
  # Manual colours
  scale_colour_manual(
    values = c(
      "0" = "blue",
      "1" = "green",
      "2" = "red",
      "3" = "black"
    )
  ) +
  
  labs(
    x = expression(italic(r)),
    y = expression(hat(r))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    panel.grid = element_blank(),
    
    # REMOVE facet labels completely
    strip.text = element_blank(),
    strip.background = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
# Save
ggsave("Output/Figure/Figure_4.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_4.png", width = 7, height = 5, bg = "white", dpi=1000)


##end Figue_4 

#using 0 to 2
#Figure_4 Recovery of r parameter##### 

# Read data
#df <- read_csv("./Data/MeanPreported_summaryNEW_NB_use2.csv")#0 to 3 Gamma(1,10)
#df <- read_csv("./Data/MeanPreported_summaryNEW_NB_Gamma.csv")# 0 to 2 Gamma(1,1)
df <- read_csv("./Data/MeanPreported_summaryNEW_NB_Gam01.csv")# 0 to 2 Gamma(0.1,0.1)

# Treat v as factor
df <- df %>%
  mutate(v = factor(v))

# Plot
ggplot(df, aes(x = p_r, y = Meanp_r, colour = v)) +
  
  # Mean estimates
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # 95% credible interval (horizontal dashes)
  geom_segment(
    aes(
      x = p_r - 0.00015, xend = p_r + 0.00015,
      y = LowerCI, yend = LowerCI
    ),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(
      x = p_r - 0.00015, xend = p_r + 0.00015,
      y = UpperCI, yend = UpperCI
    ),
    linewidth = 0.9
  ) +
  
  # True r reference lines
  geom_hline(
    aes(yintercept = p_r_true),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # Facet by ν (labels removed below)
  facet_wrap(~ v, ncol = 1) +
  
  # Manual colours
  scale_colour_manual(
    values = c(
      "0" = "blue",
      "1" = "green",
      "2" = "red"
    )
  ) +
  
  labs(
    x = expression(italic(r)),
    y = expression(hat(r))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    panel.grid = element_blank(),
    
    # REMOVE facet labels completely
    strip.text = element_blank(),
    strip.background = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
# Save
ggsave("Output/Figure/Figure_4.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_4.png", width = 7, height = 5, bg = "white", dpi=1000)




#Figure_5 Recovery of v parameter#####

# Read data
df <- read_csv("./Data/MeanV_summary3_NB_use.csv")#stops at 3

# Treat v_true as a factor
df <- df %>%
  mutate(v_true = as.factor(v_true))

# Plot
ggplot(df, aes(x = p_r, y = MeanV, color = v_true, group = v_true)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0005,
    linewidth = 0.8
  ) +
  
  # True ν reference lines
  geom_hline(
    yintercept = c(0, 1, 2, 3),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  scale_color_manual(
    values = c("#1f77b4", "#2ca02c", "#d62728", "#9467bd")#0 to 3
  ) +
  
  labs(
    x = expression(italic(r)),
    y = expression(hat(nu))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    
    panel.grid = element_blank(),  # removes all grid lines
    # remove extra space previously reserved for caption
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
# Save
ggsave("Output/Figure/Figure_5.pdf", width = 7, height = 5, bg = "white")
#end Figure_5.
ggsave("Output/Figure/Figure_5.png", width = 7, height = 5, bg = "white", dpi=1000)
#end Figure_5.


##using 0 to 2
#Figure_5 Recovery of v parameter#####

# Read data
#df <- read_csv("./Data/MeanV_summary3_NB_use2.csv")#stops at 2 Nb or norm Gammma(1,10)
#df <- read_csv("./Data/MeanV_summary3_NB_Gamma.csv")#stops at 2 Gamma(1,1)
df <- read_csv("./Data/MeanV_summary3_NB_Gam01.csv")#stops at 2 Gamma(0.1,0.1)


# Treat v_true as a factor
df <- df %>%
  mutate(v_true = as.factor(v_true))

# Plot
ggplot(df, aes(x = p_r, y = MeanV, color = v_true, group = v_true)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0005,
    linewidth = 0.8
  ) +
  
  # True ν reference lines
  geom_hline(
    yintercept = c(0, 1, 2),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  scale_color_manual(
    
    values = c("#1f77b4", "#2ca02c", "#d62728")#0 to 2
  ) +
  
  labs(
    x = expression(italic(r)),
    y = expression(hat(nu))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    
    panel.grid = element_blank(),  # removes all grid lines
    # remove extra space previously reserved for caption
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
# Save
ggsave("Output/Figure/Figure_5.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_5.png", width = 7, height = 5, bg = "white", dpi=1000)
#end Figure_5.
