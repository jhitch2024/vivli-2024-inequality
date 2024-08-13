



################################################################################
#HOUSEKEEPING
################################################################################

#loading packages
library(readr)
library(readxl)
library(pillar)
library(tidyverse)
library(GJRM)
library(wbstats)
library(summarytools)
library(tableone)

# setting working directory
setwd("...")

################################################################################
# LOADING DATA
################################################################################

# ATLAS
df_atlas_raw <- read_csv("2024_05_28 atlas_antibiotics.csv")


# GEARS
df_gears_raw <- read_excel("Venatorx surveillance data_2024_06_06.xlsx",
                       sheet = 1,
                       col_types = c("numeric", "numeric", "text", "text", "text",
                                     "text", "numeric", "text", "text", "text",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric"),
                       na = "-")

# WIID
df_wiid_raw <- read_excel("WIID_28NOV2023.xlsx")

# World Bank



  
################################################################################
# CLEANING DATA
################################################################################

df_wiid <- df_wiid_raw %>%
  filter(quality %in% c("High", "Average") & 
           areacovr == "All" & 
           popcovr == "All" & 
           resource == "Consumption")

################################################################################
# DESCRIPTIVE STATISTICS
################################################################################
# GEARS
tableone_gears <- CreateTableOne(data = df_gears)
tableone_gears$CatTable

# ATLAS
tableone_ATLAS <- CreateTableOne(data = df_atlas)
tableone_ATLAS
  

################################################################################
# DISTRIBUTIONAL COPULAS
################################################################################

eq.educ <- educ_att ~ s(age) + as.factor(hhmarstat) + as.factor(hhmale) + 
  as.factor(urban) + as.factor(num_child) + 
  as.factor(elderly) + as.factor(relig)

eq.mu <- pce.defl ~ s(age) + as.factor(hhmarstat) + as.factor(hhmale) + 
  as.factor(hheduc) + as.factor(urban) + 
  as.factor(num_child) + as.factor(elderly) + 
  as.factor(relig) + s(prov, bs = "mrf", xt = xt1, k = 15)
  
eq.si <- ~ s(age) + as.factor(hhmarstat) + 
  as.factor(num_child) + as.factor(elderly) + 
  as.factor(relig) + s(prov, bs = "mrf", xt = xt1, k = 15)

eq.theta <- ~ s(age) + as.factor(hhmarstat) + as.factor(hhmale) + 
  as.factor(hheduc) + as.factor(urban) + 
  as.factor(num_child) + as.factor(elderly) + 
  as.factor(relig) + s(prov, bs = "mrf", xt = xt1, k = 15)

form.list <- list(eq.educ, eq.mu, eq.si, eq.theta)

mod.edu <- gjrm(form.list, data = na.omit(df), ordinal = TRUE, 
                Model = "B", BivD = "N", margins = c("logit", "LN"), 
                drop.unused.levels = FALSE, gamlssfit = TRUE)

# The user first specifies the four equations for the model’s parameters of the 
# marginal distributions and of the copula, which are stored in the list 
# form.list. Continuous variables enter the model specifications via smooth 
# effects s() represented (by default) via thin-plate regression splines 
# (argument bs = "tp") with ten basis function and second order derivative 
# penalties. Spatial effects of the provinces are modeled using Markov random 
# fields with neighbourhood structure xt1 and 15 knots (argument bs = "mrf"). 
# Argument Model = "B" specifies that a bivariate model will be estimated, 
# margins = c("logit", "LN") gives the marginal distributions and BivD = "N" 
# specifies the Gaussian copula. The argument ordinal must be set to TRUE in 
# order to fit a mixed ordered-continuous model and the ordinal outcome educ_att 
# must be numeric. The optional argument gamlssfit = TRUE uses starting values 
# obtained from a univariate gamlss and drop.unused.levels = FALSE is needed 
# because not all of the provinces specified via the Markov random fields have 
# observations in the data frame df. Functions summary(), plot(), AIC() and 
# BIC() can employed in the usual manner. It is advisable to use post.check() 
# after fitting the model to produce plots of normalized quantile residuals. 
# More details, options, and the available choices for the marginal 
# distributions and copula functions can be found in the documentation of the 
# GJRM package.