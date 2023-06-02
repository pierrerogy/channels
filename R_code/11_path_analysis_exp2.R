# Path analysis
# Load libraries
library(tidyverse)
library(here)
library(brms)
library(cowplot)
library(bayestestR)
source(here::here("R_code",
                  "functions.R"))


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
                                   ## Make negative to ignore nas later in sum function
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
  ## Do it separately for number of mosquitoes
  dplyr::mutate(n_moz_log = log(n_moz + 0.01)) %>% 
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
din_treatments <-
  brms::brm(din_log ~
             exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(din_treatments)
## Check effects
bayestestR::describe_posterior(din_treatments)
## Plot
plots1a <- 
  plot_model_nice(model = din_treatments,
                  xax = "exposure",
                  yax = "din_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)

# PO4
## Fit model
po4_treatments <-
  brms::brm(po4_log ~
              exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(po4_treatments)
## Check effects
bayestestR::describe_posterior(po4_treatments)
## Plot
plots1b <- 
  plot_model_nice(model = po4_treatments,
                  xax = "exposure",
                  yax = "po4_log",
                  scale = "none",
                  type = "points",
                  data =exp2_center)


## N:P ratio
## Fit model
np_treatments <-
  brms::brm(np_log ~
              exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(np_treatments)
## Check effects
bayestestR::describe_posterior(np_treatments)
## Plot
plots1c <- 
  plot_model_nice(model = np_treatments,
                  xax = "exposure",
                  yax = "np_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
ggplot(data = model_effect,
       aes(x = exposure, 
           y = estimate__), 
       colour = subsidy) + 
  geom_point(size = 3,
             aes(colour = subsidy),
             position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = lower__, 
                    ymax = upper__,
                    colour = subsidy), 
                width = 0.2,
                position = position_dodge(0.5)) +
  geom_jitter(data = data,
              mapping = aes(x = exposure,
                            y = din_log,
                            colour = subsidy),
              alpha = 0.3) +
  ggtitle("") +
  xlab(xlab) +
  scale_x_discrete(labels = c("Exposed", "Shaded")) +
  ylab(ylab) +
  scale_color_manual(name = "Subsidy",
                     labels = c("Litter", "Litter + feces"), 
                     values = c("tan1", "tan4")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
## pH
## Fit model
pH_treatments <-
  brms::brm(pH ~
              exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.98,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(pH_treatments)
## Check effects
bayestestR::describe_posterior(pH_treatments)
## Plot
plots1d <- 
  plot_model_nice(model = pH_treatments,
                  xax = "exposure",
                  yax = "pH",
                  scale = "none",
                  type = "points",
                  data =exp2_center)
## Temperature
## Fit model
temperature_C_treatments <-
  brms::brm(temperature_C ~
              exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(temperature_C_treatments)
## Check effects
bayestestR::describe_posterior(temperature_C_treatments)
## Plot
plots1e <- 
  plot_model_nice(model = temperature_C_treatments,
                  xax = "exposure",
                  yax = "temperature_C",
                  scale = "none",
                  type = "points",
                  data =exp2_center)


# Step 2.1 - Interacting treatments on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m1 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.91,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m1)

# Bacteria
## Fit model
m101 <-
  brms::brm(bact_log  ~
              exposure*subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m101)

# Step 2.2 - Non-interacting nutrients and exposure on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m2 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m2)

# Bacteria
## Fit model
m102 <-
  brms::brm(bact_log  ~
              exposure + subsidy + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.98,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m102)

# Step 2.3 - Interacting nutrients and exposure on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m3 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(din_log_scale + po4_log_scale) + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m3)

# Bacteria
## Fit model
m103 <-
  brms::brm(bact_log  ~
              exposure*(din_log_scale + po4_log_scale) + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m103)

# Step 2.4 - Non-interacting treatments on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m4 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + din_log_scale + po4_log_scale + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m4)

# Bacteria
## Fit model
m104 <-
  brms::brm(bact_log  ~
              exposure + din_log_scale + po4_log_scale + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m104)

# Step 2.5 - Interacting treatments and nutrients on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m5 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(subsidy + din_log_scale + po4_log_scale) +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m5)

# Bacteria
## Fit model
m105 <-
  brms::brm(bact_log  ~
              exposure*(subsidy + din_log_scale + po4_log_scale)
            + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.97,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m105)

# Step 2.6 - Non-interacting nutrients and exposure on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m6 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + subsidy + din_log_scale + po4_log_scale +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.85,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m6)

# Bacteria
## Fit model
m106 <-
  brms::brm(bact_log  ~
              exposure + subsidy + din_log_scale + po4_log_scale + 
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.96,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m106)

# Step 2.7 - Interacting exposure and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m7 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*bact_log_scale + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m7)

# Bacteria
## Fit model
m107 <-
  brms::brm(bact_log  ~
              exposure*chlorophyll_ugL_log_scale + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.96,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m107)

# Step 2.8 - Non-interacting exposure and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m8 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.85,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m8)

# Bacteria
## Fit model
m108 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.89,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m108)

# Step 2.9 - Interacting nutrients and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m9 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(bact_log_scale + din_log_scale + po4_log_scale) +
                          (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m9)

# Bacteria
## Fit model
m109 <-
  brms::brm(bact_log  ~
              exposure*(chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale) +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m109)

# Step 2.10 - Non-interacting nutrients and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m10 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + din_log_scale + po4_log_scale + 
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.89,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m10)

# Bacteria
## Fit model
m110 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale + 
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m110)


# Step 2.11 - Interacting treatments and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m11 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(bact_log_scale + subsidy) +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m11)

# Bacteria
## Fit model
m111 <-
  brms::brm(bact_log  ~
              exposure*(chlorophyll_ugL_log_scale + subsidy) +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m111)

# Step 2.12 - Non-interacting treatments and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m12 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + subsidy + 
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.91,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m12)

# Bacteria
## Fit model
m112 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + subsidy + 
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m112)



# Step 2.13 - Everything interacting on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m13 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(bact_log_scale + din_log_scale + po4_log_scale + subsidy) +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m13)

# Bacteria
## Fit model
m113 <-
  brms::brm(bact_log  ~
              exposure*(chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale + subsidy) +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.91,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m113)

# Step 2.14 - Everything non-interacting  on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m14 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + din_log_scale + po4_log_scale +  subsidy +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m14)

# Bacteria
## Fit model
m114 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale + subsidy +
              (1|week) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.83,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m114)

# Step 2.15 - Picking and plotting the best model  ---------------------------------------
# Chlorophyll
## Compare fourteen models
alg_compare <- 
  brms::loo(m1, m2, m3, m4, m5,
            m6, m7, m8, m9, m10,
            m11, m12, m13, m14)
## Compute SE difference, and keep equivalent models
get_eq_models(alg_compare)
## Summary of equivalent models
summary(m11)
summary(m13)
summary(m7)
summary(m9)
summary(m12)
summary(m14)
summary(m8)
summary(m10)
## Keep only the simplest model
summary(m8)
brms::prior_summary(m8)
bayestestR::describe_posterior(m8)
## Plots
plot2a <- 
  plot_model_nice(model = m8,
                  xax = "bact_log_scale",
                  yax = "chlorophyll_ugL_log",
                  scale = "none",
                  type = "lines",
                  data =exp2_center)


# Bacteria
## Compare fourteen models
bact_compare <- 
  brms::loo(m101, m102, m103, m104, m105,
            m106, m107, m108, m109, m110,
            m111, m112, m113, m114)
## Compute SE difference, and keep equivalent models
get_eq_models(bact_compare)
## Summary of equivalent models
summary(m113)
summary(m114)
summary(m109)
## Keep only the simplest model 
summary(m114)
brms::prior_summary(m114)
bayestestR::describe_posterior(m114)
## Plots
plot2b <- 
  plot_model_nice(model = m114,
                  xax = "din_log_scale",
                  yax = "bact_log",
                  scale = "none",
                  type = "lines",
                  data =exp2_center)
plot2c <- 
  plot_model_nice(model = m114,
                  xax = "chlorophyll_ugL_log_scale",
                  yax = "bact_log",
                  scale = "none",
                  type = "lines",
                  data =exp2_center)

# Step 3  - Microorganisms and nutrients on mosquitoes ----------------------------------
# Model with interacting microorganisms and light
ma <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*(bact_log_scale*chlorophyll_ugL_log_scale) +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(ma)

# Model with non-interacting microorganisms and light
mb <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + bact_log_scale*chlorophyll_ugL_log_scale +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mb)

# Model with interacting nutrients, microorganisms and light 
mc <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*(bact_log_scale*chlorophyll_ugL_log_scale +
                          din_log_scale + po4_log_scale) +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mc)

# Model with non-interacting nutrients, microorganisms and light 
md <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + bact_log_scale*chlorophyll_ugL_log_scale +
                          din_log_scale + po4_log_scale +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(md)

# Model with interacting nutrient and light 
me <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*(din_log_scale + po4_log_scale) +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(me)

# Model with non-interacting nutrients and light 
mf <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + din_log_scale + po4_log_scale +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.85,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mf)

# Model with interacting treatments
mg <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*subsidy +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mg)

# Model with non-interacting treatments
mh <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + subsidy +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mh)

# Model with interacting light, treatment and nutrients
mi <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*(din_log_scale + po4_log_scale + subsidy) +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mi)

# Model with non-interacting light, treatment and nutrients
mj <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + din_log_scale + po4_log_scale + subsidy +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mj)

# Model with interacting microorganisms, subsidy and light
mk <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*(bact_log_scale*chlorophyll_ugL_log_scale +
                          subsidy) +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.83,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mk)

# Model with non-interacting microorganisms, subsidy and light
ml <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + bact_log_scale*chlorophyll_ugL_log_scale +
                          subsidy +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(ml)

# Model with interacting everything
mm <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure*(bact_log_scale*chlorophyll_ugL_log_scale +
                          din_log_scale + po4_log_scale + subsidy) +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mm)

# Model with non-interacting everything
mn <- 
  brms::brm(log(n_moz + 0.01)  ~
              exposure + bact_log_scale*chlorophyll_ugL_log_scale +
                          din_log_scale + po4_log_scale + subsidy +
              (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.85,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(mn)

## Compare fourteen models
mosq_compare <- 
  brms::loo(ma, mb, mc, md, me,
            mf, mg, mh, mi, mj,
            mk, ml, mm, mn)
## Compute SE difference, and keep equivalent models
get_eq_models(mosq_compare)
## Summary of equivalent models
summary(md)
summary(mm)
summary(mc)
summary(mf)
summary(mj)
summary(me)
summary(mi)
## Keep only the simplest model
summary(mf)
bayestestR::describe_posterior(mf)
brms::prior_summary(mf)
# Plot of the best model
plot3a <- 
  plot_model_nice(model = mf,
                  xax = "din_log_scale",
                  yax = "n_moz",
                  scale = "none",
                  type = "lines",
                  data =exp2_center)
plot3b <- 
  plot_model_nice(model = m114,
                  xax = "po4_log_scale",
                  yax = "n_moz",
                  scale = "none",
                  type = "lines",
                  data =exp2_center)

# Compile and save figures ------------------------------------------------
# Get legend
legend <- 
  cowplot::get_legend(plots1a)

# Figure S1
## Plot
plots1 <- 
  cowplot::plot_grid(plots1a +
                     theme(legend.position = "none") +
                     ggtitle("a"),
                   plots1b +
                     theme(legend.position = "none") +
                     ggtitle("b"),
                   plots1c +
                     theme(legend.position = "none") +
                     ggtitle("c"),
                   plots1d +
                     theme(legend.position = "none") +
                     ggtitle("d"),
                   plots1e +
                     theme(legend.position = "none") +
                     ggtitle("e"),
                   legend,
                   ncol = 2)
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figures1.jpeg"),
                plots1,
                height = 15,
                width = 10,
                bg  = "white")

# Figure 2
## Plot
plot2 <- 
  cowplot::plot_grid(plot2a +
                     theme(legend.position = "none") +
                       ggtitle("a"),
                     plot2b +
                       theme(legend.position = "none") +
                       ggtitle("b"),
                     plot2c +
                       theme(legend.position = "none") +
                       ggtitle("c"),
                     legend)
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figure2.jpeg"),
                plot2,
                height = 10,
                width = 10,
                bg  = "white")    


# Figure 3
## Plot
plot3 <- 
  cowplot::plot_grid(plot3a +
                       theme(legend.position = "none") +
                       ggtitle("a"),
                     plot3b +
                       theme(legend.position = "none") +
                       ggtitle("b"),
                     legend,
                     ncol = 3)
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figure3.jpeg"),
                plot3,
                height = 4,
                width = 11,
                bg  = "white")  
