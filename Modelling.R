
# This script is the second script of the project 
# This script is for modelling using Copula GAMLSS 

################################################################################
################################################################################

# 0. Import packages and datasets ####

# Set wd
setwd("")
dir()

# Upload packages
library(ggplot2)
library(mgcv)
library(tidyverse)
library(dplyr)
library(readxl)
library(GJRM)
library(MASS)

# Open datasets
load("cleaned_data.RData") # These are the data described at the end of script 1 "VIVLI_team1_script_cleaning_merging.R"

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# 1. Data preparation ####

# 1.1 Select variables for analysis  

df_30_full <- df_30

df_30 <- df_30_full %>%
  dplyr::select(Isolate.Id, Gender, Age.Group, Meropenem_MIC, gni_cap_usdat, j01d_beta, 
         did_total, income, Country, iso3_code, year) #Obs=7959 # Keep only variables used for analysis

# 1.2 Create datasets for different antibiotic use (complete case)

summary(is.na(df_30))

#  Beta-lactam antibiotics

df_30_beta <- df_30 %>%  
  filter(!is.na(Gender)) %>%
  filter(!is.na(gni_cap_usdat)) %>%
  filter(!is.na(j01d_beta)) #Obs=5357

summary(is.na(df_30_beta))

# Total antibiotic consumption 

df_30_did <- df_30 %>%  
  filter(!is.na(Gender)) %>%
  filter(!is.na(gni_cap_usdat)) %>%
  filter(!is.na(did_total)) #Obs=7644

summary(is.na(df_30_did))

# 1.3 Create datasets for high and low/middle income countries (with and without India)

# NOTE:
# In the final analysis we used high income countries and low/middle income countries without India.
# Therefore, the countries in the "low/middle income countries" category are actually upper and middle income countries.

# Beta-lactam antibiotics 

df_30_beta_high <- df_30_beta %>%
  filter(income == "High income")

df_30_beta_high %>%
  group_by(Country) %>%
  summarise(mean_y = mean(gni_cap_usdat)) %>%
  ungroup() %>%
  arrange(mean_y)

df_30_beta_lowmid <- df_30_beta %>%
  filter(income == "Low income" | income=="Lower middle income" | income== "Upper middle income")

df_30_beta_lowmid_noindia <- df_30_beta_lowmid %>%
  filter(Country != "India")

df_30_beta_lowmid_noindia %>%
  group_by(Country) %>%
  summarise(mean_y = mean(gni_cap_usdat)) %>%
  ungroup() %>%
  arrange(mean_y)


# Total antibiotic consumption 

df_30_did_high <- df_30_did %>%
  filter(income == "High income")

df_30_did_lowmid <- df_30_did %>%
  filter(income == "Low income" | income=="Lower middle income" | income== "Upper middle income")

df_30_did_lowmid_noindia <- df_30_did_lowmid %>%
  filter(Country != "India")

save(df_30_beta, df_30_did, file = "data_for_modelling.RData")  # These are the data used for descriptive statistics 
save(df_30_beta, df_30_beta_high, df_30_beta_lowmid, df_30_beta_lowmid_noindia,
     df_30_did, df_30_did_high, df_30_did_lowmid, df_30_beta_lowmid_noindia,file = "full_data_for_modelling.RData")


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# 2. Run Copula GAMLSS   ####

load("full_data_for_modelling.RData")

# 2.1 Model specification 

# 2.1.1 Beta-lactam antibiotics 

# (1) High income

eq.amr <- Meropenem_MIC ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)  
eq.au <- j01d_beta ~ Age.Group + Gender  + year + s(gni_cap_usdat, bs="ps", k=4)
eq.si <- ~s(gni_cap_usdat, bs="ps", k=4)
eq.theta <- ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)
form.list <- list(eq.amr, eq.au, eq.si, eq.theta)

m_beta_high <- gjrm(form.list, data = df_30_beta_high, model = "B", margins = c("logit", "LN"),
           drop.unused.levels = TRUE)


# (2) Low and middle income excluding India
eq.amr <- Meropenem_MIC ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)  
eq.au <- j01d_beta ~ Age.Group + Gender  + year + s(gni_cap_usdat, bs="ps", k=4)
eq.si <- ~s(gni_cap_usdat, bs="ps", k=4)
eq.theta <- ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)
form.list <- list(eq.amr, eq.au, eq.si, eq.theta)

m_beta_lowmid <- gjrm(form.list, data = df_30_beta_lowmid_noindia, model = "B", margins = c("logit", "LN"),
           drop.unused.levels = TRUE)

# 2.1.2 Total antibiotic consumption 

# (1) High income

eq.amr <- Meropenem_MIC ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)  
eq.au <- did_total ~ Age.Group + Gender  + year + s(gni_cap_usdat, bs="ps", k=4)
eq.si <- ~s(gni_cap_usdat, bs="ps", k=4)
eq.theta <- ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)
form.list <- list(eq.amr, eq.au, eq.si, eq.theta)

m_did_high <- gjrm(form.list, data = df_30_did_high, model = "B", margins = c("logit", "LN"),
               drop.unused.levels = TRUE)

#(2) Low and middle income countries excluding India 

eq.amr <- Meropenem_MIC ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)  
eq.au <- did_total ~ Age.Group + Gender  + year + s(gni_cap_usdat, bs="ps", k=4)
eq.si <- ~s(gni_cap_usdat, bs="ps", k=4)
eq.theta <- ~ Age.Group + Gender + year + s(gni_cap_usdat, bs="ps", k=4)
form.list <- list(eq.amr, eq.au, eq.si, eq.theta)

m_did_lowmid <- gjrm(form.list, data = df_30_did_lowmid_noindia, model = "B", margins = c("logit", "LN"),
           drop.unused.levels = TRUE)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# 2.2 Plotting using MIC & Beta-lactam model ####

# 2.2.1 High-income countries 

# (1) Plot for MIC

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_beta_high$gni_cap_usdat), max(df_30_beta_high$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_beta_high$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_beta_high$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_beta_high$year)

predictions <- predict(m_beta_high$gam1, newdata = newdata, type = "terms", se.fit = TRUE)

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)


ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "High-Income Countries",
       x = "Gross National Income per Capita (Current USD)",
       y = "Partial Effect of GNI on MIC") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
  scale_x_continuous(breaks = c(20000, 40000, 60000)) 


# (2) Plot for Beta-lactam 

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_beta_high$gni_cap_usdat), max(df_30_beta_high$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_beta_high$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_beta_high$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_beta_high$year)

predictions <- predict(m_beta_high$gam2, newdata = newdata, type = "terms", se.fit = TRUE)

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)


# Deal with 0 in y-aix shown as 0e+00
library(scales)

custom_format <- function(x) {
  # Replace 0 with "0" and format other numbers with 2 decimal places
  ifelse(x == 0, "0", format(x, scientific = FALSE, digits = 2))
}

ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "High-Income Countries",
       x = "Gross National Income per Capita (Current USD)",
       y = "Partial Effect of GNI on Consumption of other beta-lactams") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(breaks = c(20000, 40000, 60000)) +
  scale_y_continuous(labels = custom_format)


# 2.2.2 Low and middle income countries (excluding India)

# (1) Plot for MIC

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_beta_lowmid_noindia$gni_cap_usdat), max(df_30_beta_lowmid_noindia$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_beta_lowmid_noindia$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_beta_lowmid_noindia$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_beta_lowmid_noindia$year)

predictions <- predict(m_beta_lowmid$gam1, newdata = newdata, type = "terms", se.fit = TRUE)

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)


ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "Upper Middle-Income Countries",
       x = "Gross National Income per Capita (Current USD)",
       y = "Partial Effect of GNI on MIC") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  scale_y_continuous(labels = custom_format)
# +
#   scale_x_continuous(breaks = c(2000, 6000, 10000)) 


# (2) Plot for Beta-lactam 

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_beta_lowmid_noindia$gni_cap_usdat), max(df_30_beta_lowmid_noindia$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_beta_lowmid_noindia$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_beta_lowmid_noindia$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_beta_lowmid_noindia$year)

predictions <- predict(m_beta_lowmid$gam2, newdata = newdata, type = "terms", se.fit = TRUE)

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

custom_format <- function(x) {
  # Replace 0 with "0" and format other numbers with 2 decimal places
  ifelse(x == 0, "0", format(x, scientific = FALSE, digits = 2))
}

ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "Upper Middle-Income Countries",
       x = "Gross National Income per Capita (Current USD)",
       y = "Partial Effect of GNI on Beta-Lactam Antibiotic Consumption") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  #scale_x_continuous(breaks = c(2000, 6000, 10000)) +
  scale_y_continuous(labels = custom_format)

plot(df_30_beta_lowmid_noindia$gni_cap_usdat, df_30_beta_lowmid_noindia$did_total)


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# 2.3 Plotting using MIC & total antibiotic consumption model ####

# 2.3.1 High-income countries 

# (1) Plot for MIC

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_did_high$gni_cap_usdat), max(df_30_did_high$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_did_high$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_did_high$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_did_high$year)

predictions <- predict(m_did_high$gam1, newdata = newdata, type = "terms", se.fit = TRUE)

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)


ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "High-Income Countries",
       x = "Gross National Income per Capita (USD)",
       y = "Partial Effect of GNI on MIC") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  scale_x_continuous(breaks = c(20000, 40000, 60000)) 



# (2) Plot for total antibiotic consumption 

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_did_high$gni_cap_usdat), max(df_30_did_high$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_did_high$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_did_high$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_did_high$year)

predictions <- predict(m_did_high$gam2, newdata = newdata, type = "terms", se.fit = TRUE)

predictions_response <- predict(m_did_high$gam2, newdata = newdata, type = "response", se.fit = TRUE)
predictions_response

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)


ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "High-Income Countries",
       x = "Gross National Income per Capita (USD)",
       y = " Partial Effect of GNI on Total Antibiotic Consumption") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(breaks = c(20000, 40000, 60000)) +
  scale_y_continuous(labels = custom_format)





# 2.3.2 Low and middle income countries

# (1) Plot for MIC

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_did_lowmid_noindia$gni_cap_usdat), max(df_30_did_lowmid_noindia$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_did_lowmid_noindia$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_did_lowmid_noindia$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_did_lowmid_noindia$year)

predictions <- predict(m_did_lowmid$gam1, newdata = newdata, type = "terms", se.fit = TRUE)

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)


ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "Upper Middle-Income Countries",
       x = "Gross National Income per Capita (USD)",
       y = "Partial Effect of GNI on MIC") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)
  ) 
# +
#   scale_x_continuous(breaks = c(2000, 6000, 10000)) 




# (2) Plot for total antibiotic consumption 

newdata <- data.frame(gni_cap_usdat = seq(min(df_30_did_lowmid_noindia$gni_cap_usdat), max(df_30_did_lowmid_noindia$gni_cap_usdat), length.out = 100))

newdata$Age.Group <- df_30_did_lowmid_noindia$Age.Group[1]  # Assigning a representative value or level
newdata$Gender <- df_30_did_lowmid_noindia$Gender[1]        # Assigning a representative value or level
newdata$year <- mean(df_30_did_lowmid_noindia$year)

predictions <- predict(m_did_lowmid$gam2, newdata = newdata, type = "terms", se.fit = TRUE)
predictions$fit

predictions_response <- predict(m_did_lowmid$gam2, newdata = newdata, type = "response", se.fit = TRUE)
predictions_response

fitted_values <- predictions$fit[, "s(gni_cap_usdat)"]
se_fit <- predictions$se.fit[, "s(gni_cap_usdat)"]

lower_bound <- fitted_values - 1.96 * se_fit
upper_bound <- fitted_values + 1.96 * se_fit

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

plot_data <- data.frame(
  gni_cap_usdat = newdata$gni_cap_usdat,
  fit = fitted_values,
  lower = lower_bound,
  upper = upper_bound
)

ggplot(plot_data, aes(x = gni_cap_usdat, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "lightblue") +
  labs(title = "Upper Middle-Income Countries",
       x = "Gross National Income per Capita (USD)",
       y = "Predicted Total Antibiotic Consumption") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.text.x = element_text(angle = 0, hjust = 1), # Optional: angle x-axis labels if needed
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  #scale_x_continuous(breaks = c(2000, 6000, 10000)) +
  scale_y_continuous(labels = custom_format)


################################################################################
# Plot of total antibiotic consumption over time grouped by income group
df <- rbind(df_30_did_high, df_30_did_lowmid_noindia)

df_plot <- df %>%
  group_by(income, year) %>%
  summarise(mean_did = mean(did_total)) %>%
  ungroup()

p <- ggplot(df_plot, aes(year, mean_did, color = income)) +
  labs(color = "Income Classification",
       title = "Total antibiotic consumption by income group",
       x = "Year",
       y = "Total antibiotic consumption (DDD/1,000/day)") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_y_continuous(labels = custom_format) +
  scale_x_continuous(breaks = seq(2006, 2019, by = 1)) +
  geom_line(size = 1.25)
p

df <- rbind(df_30_beta_high, df_30_beta_lowmid_noindia)

df_plot <- df %>%
  group_by(income, year) %>%
  summarise(n = n(),
            mean_j01d = mean(j01d_beta)) %>%
  ungroup()

ggplot(df_plot, aes(year, mean_j01d, color = income)) +
  labs(color = "Income Classification",
       title = "J01D Other beta-lactams consumption by income group",
       x = "Year",
       y = "J01D Other lactams consumption (DDD/1,000/day)") +
  theme_minimal(base_size = 15) +  # Increase base_size for better readability
  theme(panel.background = element_rect(fill = "white"), # Set background to white
        panel.grid.major = element_blank(),              # Remove major grid lines
        panel.grid.minor = element_blank(),              # Remove minor grid lines
        axis.ticks.x = element_line(),                    # Optionally keep x-axis ticks
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15)),  # Adjust space between x-axis title and axis
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_y_continuous(labels = custom_format) +
  scale_x_continuous(breaks = seq(2006, 2019, by = 1)) +
  geom_line(size = 1.25)