# Path analysis
# Load libraries
library(tidyverse)
library(here)
library(lme4)
library(MCMCglmm)
library(ggeffects)
library(DHARMa)
library(afex)
library(ggplot2)
source(here::here("R_code",
                  "functions.R"))
library(mgcv)


# Prepare and transform data ----------------------------------------------
# Load mosquito data
mozzies <- 
  read.csv(here::here("microchannels",
                      "appdata", 
                      "mosquitoes.csv")) %>% 
  ## Add average of two wings
  dplyr::mutate(wing_length = (left_wing_mm + right_wing_mm)/2)

# Load chemical data
exp2 <- 
  read.csv(here::here("microchannels",
                      "appdata", 
                      "weekly_measurements_exp2.csv")) %>% 
  ## Make factors for models
  dplyr::mutate(across(c(larvae, subsidy, shading, 
                         cup_number, week),
                       ~as.factor(.))) %>% 
  ## Add DIN column
  dplyr::mutate(din = nh4 +  no2 + no3) %>% 
  ## Column with number of mosquitoes, 0 for absent and first week of present
  dplyr::mutate(n_moz = ifelse(larvae == "absent",
                               0, ifelse(week == 0,
                                         0,  6))) %>% 
  ## Get cumulative number of dead or pupating mosquitoes each week
  dplyr::left_join(mozzies %>% 
                     dplyr::select(cup_number, week_of_departure) %>% 
                     dplyr::group_by(cup_number, week_of_departure) %>% 
                     dplyr::tally() %>% 
                     dplyr::mutate(cumdeath = cumsum(n),
                                   ##ake negative to ignore nas later in sum function
                                   cumdeath = -cumdeath) %>% 
                     dplyr::ungroup() %>% 
                     ## if a mosquito dies a given week, only the next week measurement
                     ## will have one less mosquito
                     dplyr::mutate(week = week_of_departure + 1,
                                   week = as.factor(week),
                                   cup_number = as.factor(cup_number)) %>% 
                     dplyr::select(-n, - week_of_departure),
                   by = c("cup_number", "week")) %>% 
  ## Now compute how many mosquitoes in a cup each week
  dplyr::rowwise() %>% 
  dplyr::mutate(n_moz = sum(c(n_moz, cumdeath), 
                            na.rm = TRUE)) %>% 
  dplyr::select(-cumdeath) %>% 
  ## Lets not forget the misplaced one
  dplyr::mutate(n_moz = ifelse(cup_number == 46 & week %in% c("1", "2"),
                               1, ifelse(cup_number == 46 & week == "3",
                                         0, n_moz))) %>% 
  ## Add N:P ratio
  dplyr::mutate(np = din/po4) %>% 
  ## Ungroup
  dplyr::ungroup() %>% 
  ## Put variables on log scale
  dplyr::mutate(across(c(din, po4, np,
                         chlorophyll_ugL, bact),
                       ~ log(.x),
                       .names = "{.col}_log"))%>% 
  ## Fix one T value that was mistyped
  dplyr::mutate(temperature_C = ifelse(temperature_C == 13.9,
                                       23.9, temperature_C)) %>% 
  ## Remove all NAs
  dplyr::filter(!is.na(din) & !is.na(po4) &
                  !is.na(bact) & !is.na(chlorophyll_ugL))
# Scale data 
exp2_center <- 
  exp2 %>% 
  ## Rename shading column for next action (contains "din")
  dplyr::rename(exposure = shading) %>% 
  ## Center around 0
  dplyr::mutate(across(contains(c("din","po4","np","chlorophyll_ugL","bact")),
                       ~as.numeric(scale(.)),
                       .names = "{.col}_scale")) 

# Step 1 - Effects of nonmoz treatments on nutrients, pH and T-----------------------------
# DIN
## Fit model
din_exp2_reslight <-
  MCMCglmm::MCMCglmm(din_log ~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian",                      
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(din_exp2_reslight$Sol)
autocorr.diag(din_exp2_reslight$VCV)
## Test
summary(din_exp2_reslight)
## Plot
plot1 <- 
  plot_model_nice(model = din_exp2_reslight,
                  parameter = "din",
                  scale = "log",
                  type = "points")

# PO4
## Fit model
po4_exp2_reslight <-
  MCMCglmm::MCMCglmm(po4_log ~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian",                      
                     nitt = 100000, 
                     burnin = 10000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(po4_exp2_reslight$Sol)
autocorr.diag(po4_exp2_reslight$VCV)
## Test
summary(po4_exp2_reslight)
## Plot
plot2 <- 
  plot_model_nice(model = po4_exp2_reslight ,
                  parameter = "po4",
                  scale = "log",
                  type = "points")

## N:P ratio
## Fit model
np_exp2_reslight <-
  MCMCglmm::MCMCglmm(np_log ~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian",                      
                     nitt = 100000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(np_exp2_reslight$Sol)
autocorr.diag(np_exp2_reslight$VCV)
## Test
summary(np_exp2_reslight)
## Plot
plot3 <- 
  plot_model_nice(model = np_exp2_reslight,
                  parameter = "np",
                  scale = "log",
                  type = "points")

## pH
## Fit model
pH_exp2_reslight <-
  MCMCglmm::MCMCglmm(pH~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(pH_exp2_reslight$Sol)
autocorr.diag(pH_exp2_reslight$VCV)
## Test
summary(pH_exp2_reslight)
## Plot
plots1 <- 
  plot_model_nice(model = pH_exp2_reslight,
                  parameter = "pH",
                  type = "points")

## Temperature
## Fit model
temperature_exp2_reslight <-
  MCMCglmm::MCMCglmm(temperature_C~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 100000, 
                     burnin = 1000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(temperature_exp2_reslight$Sol)
autocorr.diag(temperature_exp2_reslight$VCV)
## Test
summary(temperature_exp2_reslight)
## Plot
plots2 <- 
  plot_model_nice(model = temperature_exp2_reslight ,
                  parameter = "temperature_C",
                  type = "points")


# Step 2.1 - Treatments on microorganisms ------------------------------------
# Clhorophyll
## Fit model
chloro_exp2_treatlight <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 100000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_treatlight$Sol)
autocorr.diag(chloro_exp2_treatlight$VCV)

# Bacteria
## Fit model
bact_exp2_treatlight <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 5000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_treatlight$Sol)
autocorr.diag(bact_exp2_treatlight$VCV)

# Step 2.2 - Nutrients on microorganisms ------------------------------------
# Chlorophyll
## Fit model
chloro_exp2_nutlight <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*(din_log_scale + po4_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 100000, 
                     burnin = 1000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_nutlight$Sol)
autocorr.diag(chloro_exp2_nutlight$VCV)

# Bacteria 
## Fit model
bact_exp2_nutlight <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*(din_log_scale + po4_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 10000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_nutlight$Sol)
autocorr.diag(bact_exp2_nutlight$VCV)

# Step 2.3 - Treatments and nutrients on microorganisms ------------------------------------
# Chlorophyll
## Fit model
chloro_exp2_treatnutlight <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*(subsidy + din_log_scale + po4_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 120000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_treatnutlight$Sol)
autocorr.diag(chloro_exp2_treatnutlight$VCV)

# Bacteria
## Fit model
bact_exp2_treatnutlight <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*(subsidy + din_log_scale + po4_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 5000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_treatnutlight$Sol)
autocorr.diag(bact_exp2_treatnutlight$VCV)

# Step 2.4 - Intra-microorganisms -------------------------------------------
## Fit model
chloro_exp2_bact <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*bact_log_scale,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 100000, 
                     burnin = 10000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_bact$Sol)
autocorr.diag(chloro_exp2_bact$VCV)

## Fit model
bact_exp2_chloro <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*chlorophyll_ugL_log_scale,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 100000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_chloro$Sol)
autocorr.diag(bact_exp2_chloro$VCV)

# Step 2.5 - Intra-microorganisms and nutrients -------------------------------------------
# Chlorophyll
## Fit model
chloro_exp2_bactnut <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*(din_log_scale + po4_log_scale + 
                                bact_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 2000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_bactnut$Sol)
autocorr.diag(chloro_exp2_bactnut$VCV)

# Bacteria
## Fit model
bact_exp2_chloronut <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*(din_log_scale + po4_log_scale +
                                  chlorophyll_ugL_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 1000,
                     thin = 200,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_chloronut$Sol)
autocorr.diag(bact_exp2_chloronut$VCV)

# Step 2.6 - Intra-microorganisms and treatment -------------------------------------------
# Chlorophyll
## Fit model
chloro_exp2_bacttreat <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*(subsidy + bact_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 100000, 
                     burnin = 15000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_bacttreat$Sol)
autocorr.diag(chloro_exp2_bacttreat$VCV)

# Bacteria
## Fit model
bact_exp2_chlorotreat <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*(subsidy + chlorophyll_ugL_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 1000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_chlorotreat$Sol)
autocorr.diag(bact_exp2_chlorotreat$VCV)

# Step 2.7 - Everything on microorganisms ------------------------------------
# Chorophyll
## Fit model
chloro_exp2_all <-
  MCMCglmm::MCMCglmm(chlorophyll_ugL_log ~
                       exposure*(subsidy + din_log_scale +
                                  po4_log_scale + bact_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 120000, 
                     burnin = 10000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_all$Sol)
autocorr.diag(chloro_exp2_all$VCV)

# Bacteria
## Fit model
bact_exp2_all <-
  MCMCglmm::MCMCglmm(bact_log ~
                       exposure*(subsidy + din_log_scale +
                                  po4_log_scale + chlorophyll_ugL_log_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 15000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_all$Sol)
autocorr.diag(bact_exp2_all$VCV)

# Step 2.8 - Picking and plotting the best model ---------------------------------------
# Chlorophyll
## Highest to lowest
summary(chloro_exp2_nutlight)$DIC
summary(chloro_exp2_treatlight)$DIC
summary(chloro_exp2_treatnutlight)$DIC
summary(chloro_exp2_bact)$DIC
summary(chloro_exp2_bactnut)$DIC
summary(chloro_exp2_all)$DIC
summary(chloro_exp2_bacttreat)$DIC
## Summary
summary(chloro_exp2_bacttreat)
## Plots
plot45 <- 
  plot_model_nice(model = chloro_exp2_bacttreat,
                  parameter = "bact_log_scale",
                  scale = "log",
                  type = "lines")


## Bacteria
## Highest to lowest
summary(bact_exp2_chloro)$DIC
summary(bact_exp2_treatlight)$DIC
summary(bact_exp2_nutlight)$DIC
summary(bact_exp2_chlorotreat)$DIC
summary(bact_exp2_chloronut)$DIC
summary(bact_exp2_treatnutlight)$DIC
summary(bact_exp2_all)$DIC
## Summary
summary(bact_exp2_all)
## Plot
## DIN
plot67 <- 
  plot_model_nice(model = bact_exp2_all,
                  parameter = "din_log_scale",
                  scale = "log",
                  type = "lines")
plot89 <- 
  plot_model_nice(model = bact_exp2_all,
                  parameter = "po4_log_scale",
                  scale = "log",
                  type = "lines")
plot1011 <- 
  plot_model_nice(model = bact_exp2_all,
                  parameter = "chlorophyll_ugL_log_scale",
                  scale = "log",
                  type = "lines")

# Step 3  - Microorganisms and nutrients on mosquitoes ----------------------------------
# Model with microorganisms and light
mod_gam1 <- 
  mgcv::gam(n_moz ~ s(bact_log_scale, k = 15) + s(chlorophyll_ugL_log_scale, k = 15) + 
              exposure + s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15) + 
              s(bact_log_scale, k = 15, by = exposure) + s(chlorophyll_ugL_log_scale, k = 15, by = exposure) + 
              s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15, by = exposure) + 
              s(cup_number, bs = "re") + s(week, bs = "re"),
            data = exp2_center, method="REML")

# Model with nutrientsand light
mod_gam2 <- 
  mgcv::gam(n_moz ~ s(din_log_scale, k = 15) + s(po4_log_scale, k = 15) +
              s(din_log_scale, k = 15, by = exposure) + 
              s(po4_log_scale, k = 15, by = exposure) + exposure + 
              s(cup_number, bs = "re") + s(week, bs = "re"), 
            data = exp2_center, method="REML")

# Model with treatment and light
mod_gam3 <- 
  mgcv::gam(n_moz ~ subsidy*exposure +
               s(cup_number, bs = "re") + s(week, bs = "re"), 
             data = exp2_center, method="REML")


# Model with microorganisms and light and nutrients
mod_gam4 <- 
  mgcv::gam(n_moz ~ exposure + s(din_log_scale, k = 15) + s(po4_log_scale, k = 15) +
              s(din_log_scale, k = 15, by = exposure) + s(po4_log_scale, k = 15, by = exposure) + 
              s(bact_log_scale, k = 15) + s(chlorophyll_ugL_log_scale, k = 15) + 
              s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15) + 
              s(bact_log_scale, k = 15, by = exposure) + s(chlorophyll_ugL_log_scale, k = 15, by = exposure)+ 
              s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15, by = exposure) + 
              s(cup_number, bs = "re") + s(week, bs = "re"), 
            data = exp2_center, 
            method="REML")

# Model with treatment and light and nutrient
mod_gam5 <- 
  mgcv::gam(n_moz ~ s(din_log_scale, k = 15) + s(po4_log_scale, k = 15) +
              s(din_log_scale, k = 15, by = exposure) + 
              s(po4_log_scale, k = 15, by = exposure) + subsidy*exposure + 
              s(cup_number, bs = "re") + s(week, bs = "re"), 
            data = exp2_center, 
            method="REML")

# Model with microorganisms and light and treatment
mod_gam6 <- 
  mgcv::gam(n_moz ~ s(bact_log_scale, k = 15) + s(chlorophyll_ugL_log_scale, k = 15) + 
              exposure*subsidy + s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15) + 
              s(bact_log_scale, k = 15, by = exposure) + s(chlorophyll_ugL_log_scale, k = 15, by = exposure)+ 
              s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15, by = exposure) + 
              s(cup_number, bs = "re") + s(week, bs = "re"), 
            data = exp2_center, 
            method="REML")

# Model with microorganisms and light and nutrients and treatment
mod_gam7 <- 
  mgcv::gam(n_moz ~ s(din_log_scale, k = 15) + s(po4_log_scale, k = 15) + 
         s(din_log_scale, k = 15, by = exposure) + s(po4_log_scale, k = 15, by = exposure) + 
         s(bact_log_scale, k = 15) + s(chlorophyll_ugL_log_scale, k = 15) + 
         exposure*subsidy + s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15) + 
         s(bact_log_scale, k = 15, by = exposure) + s(chlorophyll_ugL_log_scale, k = 15, by = exposure)+ 
         s(bact_log_scale, chlorophyll_ugL_log_scale, k = 15, by = exposure) + 
         s(cup_number, bs = "re") + s(week, bs = "re"), 
      data = exp2_center, 
      method="REML")

# Picking the best model
## Ranked from highest to lowest
AIC(mod_gam2)
AIC(mod_gam5)
AIC(mod_gam6)
AIC(mod_gam1)
AIC(mod_gam7)
AIC(mod_gam4)
AIC(mod_gam3)



# Summary of the best model
summary(mod_gam3)
anova(mod_gam3)

# Plot of the best model
plot_model_nice(model = mod_gam3,
                parameter = "subsidy",
                scale = "none",
                type = "points")
        