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
  dplyr::mutate(across(contains(c("din","po4","np","chlorophyll_ugL","bact", "n_moz")),
                       ~as.numeric(scale(.)),
                       .names = "{.col}_scale"))
  

# Step 0 - Effects of treatments on pH and T-----------------------------
## pH
## Fit model
pH_treatments <-
  brms::brm(pH ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.97,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(pH_treatments)
bayestestR::describe_posterior(pH_treatments)
## Check effects
all_tests(model = pH_treatments)
## Plot
plots1a <- 
  plot_model_nice(model = pH_treatments,
                  xax = "exposure",
                  yax = "pH",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
plots1b <- 
  plot_model_nice(model = pH_treatments,
                  xax = "n_moz_log_scale",
                  yax = "pH",
                  scale = "none",
                  type = "combo",
                  data = exp2_center)

## Temperature
## Fit model
temperature_C_treatments <-
  brms::brm(temperature_C ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(temperature_C_treatments)
bayestestR::describe_posterior(temperature_C_treatments)
## Check effects
all_tests(model = temperature_C_treatments)
## Plot
plots1c <- 
  plot_model_nice(model = temperature_C_treatments,
                  xax = "exposure",
                  yax = "temperature_C",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
plots1d <- 
  plot_model_nice(model = temperature_C_treatments,
                  xax = "n_moz_log_scale",
                  yax = "temperature_C",
                  scale = "none",
                  type = "combo",
                  data = exp2_center)

# Step 1.1  - Treatments, microorganisms and mosquitos on nutrients ----------------------------------
# DIN
## Fit model
din_treatments <-
  brms::brm(din_log ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(din_treatments)
bayestestR::describe_posterior(din_treatments)
## Test
all_tests(model = din_treatments)
## Plot
plot2a <- 
  plot_model_nice(model = din_treatments,
                  xax = "exposure",
                  yax = "din_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
plot2b <- 
  plot_model_nice(model = din_treatments,
                  xax = "n_moz_log_scale",
                  yax = "din_log",
                  scale = "none",
                  type = "combo",
                  data = exp2_center)

# PO4
## Fit model
po4_treatments <-
  brms::brm(po4_log ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
# Check assumptions
plot(po4_treatments)
bayestestR::describe_posterior(po4_treatments)
## Test
all_tests(model = po4_treatments)
## Plot
plot2c <- 
  plot_model_nice(model = po4_treatments,
                  xax = "exposure",
                  yax = "po4_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)

# N:P ratio
## Fit model
np_treatments <-
  brms::brm(np_log ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(np_treatments)
bayestestR::describe_posterior(np_treatments)
## Test
all_tests(model = np_treatments)
## Plots
plot2d <- 
  plot_model_nice(model = np_treatments,
                  xax = "exposure",
                  yax = "np_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
plot2e <- 
  plot_model_nice(model = np_treatments,
                  xax = "n_moz_log_scale",
                  yax = "np_log",
                  scale = "none",
                  type = "combo",
                  data = exp2_center)

# Step 2.1 - Interacting treatments on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m1 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.91,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m1)
bayestestR::describe_posterior(m1)

# Bacteria
## Fit model
m101 <-
  brms::brm(bact_log  ~
              exposure*subsidy*n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m101)
bayestestR::describe_posterior(m101)

# Step 2.2 - Non-interacting nutrients and exposure on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m2 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + subsidy + n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m2)
bayestestR::describe_posterior(m2)

# Bacteria
## Fit model
m102 <-
  brms::brm(bact_log  ~
              exposure + subsidy + n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.98,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m102)
bayestestR::describe_posterior(m102)

# Step 2.3 - Interacting nutrients and exposure on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m3 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(din_log_scale + po4_log_scale + n_moz_log_scale) + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m3)
bayestestR::describe_posterior(m3)

# Bacteria
## Fit model
m103 <-
  brms::brm(bact_log  ~
              exposure*(din_log_scale + po4_log_scale + n_moz_log_scale) + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m103)
bayestestR::describe_posterior(m103)

# Step 2.4 - Non-interacting treatments on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m4 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + din_log_scale + po4_log_scale + n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m4)
bayestestR::describe_posterior(m4)

# Bacteria
## Fit model
m104 <-
  brms::brm(bact_log  ~
              exposure + din_log_scale + po4_log_scale + n_moz_log_scale + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m104)
bayestestR::describe_posterior(m104)

# Step 2.5 - Interacting treatments and nutrients on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m5 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(subsidy + din_log_scale + po4_log_scale + n_moz_log_scale) +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m5)
bayestestR::describe_posterior(m5)

# Bacteria
## Fit model
m105 <-
  brms::brm(bact_log  ~
              exposure*(subsidy + din_log_scale + po4_log_scale + n_moz_log_scale)
            + (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.97,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m105)
bayestestR::describe_posterior(m105)

# Step 2.6 - Non-interacting nutrients and exposure on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m6 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + subsidy + din_log_scale + po4_log_scale + n_moz_log_scale +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.85,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m6)
bayestestR::describe_posterior(m6)

# Bacteria
## Fit model
m106 <-
  brms::brm(bact_log  ~
              exposure + subsidy + din_log_scale + po4_log_scale + n_moz_log_scale + 
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.96,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m106)
bayestestR::describe_posterior(m106)

# Step 2.7 - Interacting nutrients and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m7 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(bact_log_scale + din_log_scale + n_moz_log_scale + po4_log_scale + n_moz_log_scale) +
                          (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m7)
bayestestR::describe_posterior(m7)

# Bacteria
## Fit model
m107 <-
  brms::brm(bact_log  ~
              exposure*(chlorophyll_ugL_log_scale +  din_log_scale + po4_log_scale+ n_moz_log_scale) +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m107)
bayestestR::describe_posterior(m107)

# Step 2.8 - Non-interacting nutrients and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m8 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + din_log_scale + po4_log_scale + n_moz_log_scale + 
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.89,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m8)
bayestestR::describe_posterior(m8)

# Bacteria
## Fit model
m108 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale + n_moz_log_scale + 
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m108)
bayestestR::describe_posterior(m108)


# Step 2.9 - Interacting treatments and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m9 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(bact_log_scale + subsidy + n_moz_log_scale) +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m9)
bayestestR::describe_posterior(m9)

# Bacteria
## Fit model
m109 <-
  brms::brm(bact_log  ~
              exposure*(chlorophyll_ugL_log_scale + subsidy + n_moz_log_scale) +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m109)
bayestestR::describe_posterior(m109)

# Step 2.10 - Non-interacting treatments and microorganisms on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m10 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + subsidy + n_moz_log_scale + 
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.91,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m10)
bayestestR::describe_posterior(m10)

# Bacteria
## Fit model
m110 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + subsidy + n_moz_log_scale + 
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.93,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m110)
bayestestR::describe_posterior(m110)

# Step 2.11 - Everything interacting on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m11 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure*(bact_log_scale + din_log_scale + po4_log_scale + subsidy + n_moz_log_scale) +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m11)
bayestestR::describe_posterior(m11)

# Bacteria
## Fit model
m111 <-
  brms::brm(bact_log  ~
              exposure*(chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale + subsidy + n_moz_log_scale) +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.91,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m111)
bayestestR::describe_posterior(m111)

# Step 2.12 - Everything non-interacting  on microorganisms ------------------------------------
# Chlorophyll
## Fit model
m12 <-
  brms::brm(chlorophyll_ugL_log  ~
              exposure + bact_log_scale + din_log_scale + po4_log_scale +  subsidy + n_moz_log_scale +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.87,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m12)
bayestestR::describe_posterior(m12)

# Bacteria
## Fit model
m112 <-
  brms::brm(bact_log  ~
              exposure + chlorophyll_ugL_log_scale + din_log_scale + po4_log_scale + subsidy + n_moz_log_scale +
              (1|week) + (1|cup_number),
            iter = 5000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.85,
                           max_treedepth = 10),
            data = exp2_center)
## Check assumptions
plot(m112)
bayestestR::describe_posterior(m112)

# Step 2.13 - Picking and plotting the best model  ---------------------------------------
# Chlorophyll
## Compare twelve models
alg_compare <- 
  brms::loo(m1, m2, m3, m4, m5,
            m6, m7, m8, m9, m10,
            m11, m12)
## Compute SE difference, and keep equivalent models
get_eq_models(alg_compare)
## Summary of equivalent models
summary(m10)
summary(m1)
summary(m8)
summary(m12)
summary(m9)
summary(m7)
summary(m11)
## Keep only the simplest model
summary(m1)
all_tests(model = m1)
## Plots
plot3a <- 
  plot_model_nice(model = m1,
                  xax = "exposure",
                  yax = "chlorophyll_ugL_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
plot3b <- 
  plot_model_nice(model = m1,
                  xax = "n_moz_log_scale",
                  yax = "chlorophyll_ugL_log",
                  scale = "none",
                  type = "combo",
                  data = exp2_center)

# Bacteria
## Compare fourteen models
bact_compare <- 
  brms::loo(m101, m102, m103, m104, m105,
            m106, m109, m110, m109, m110, 
            m111, m112)
## Compute SE difference, and keep equivalent models
get_eq_models(bact_compare)
## Summary of equivalent models
summary(m112)
summary(m111)
summary(m105)
## Keep only the simplest model 
summary(m105)
all_tests(model = m105)
## Plots
plot3c <- 
  plot_model_nice(model = m105,
                  xax = "exposure",
                  yax = "bact_log",
                  scale = "none",
                  type = "points",
                  data = exp2_center)
plot3d <- 
  plot_model_nice(model = m105,
                  xax = "din_log_scale",
                  yax = "bact_log",
                  scale = "none",
                  type = "din_log_scale",
                  data = exp2_center)

# Compile and save figures ------------------------------------------------
# Get legend
legend <- 
  cowplot::get_legend(plots1a)

# Figure S1
## Plot
plots1 <- 
  cowplot::plot_grid(plots1a +
                       theme(legend.position = "none") +
                       ggtitle("(a)"),
                     plots1b +
                       theme(legend.position = "none") +
                       ggtitle("(b)"),
                     plots1c +
                       theme(legend.position = "none") +
                       ggtitle("(c)"),
                     plots1d +
                       theme(legend.position = "none") +
                       ggtitle("(d)"),
                     ncol = 2)
## Combine with legend
plots1 <- 
  cowplot::plot_grid(legend,
                     plots1,
                     nrow = 2,
                     rel_heights = c (0.1, 0.9)) 
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figures1.jpeg"),
                plots1,
                height = 11,
                width = 10,
                bg  = "white")

# Figure 2
## Plot
plot2 <- 
  cowplot::plot_grid(plot2a +
                     theme(legend.position = "none") +
                       ggtitle("(a)"),
                     plot2b +
                       theme(legend.position = "none") +
                       ggtitle("(b)"),
                     plot2c +
                       theme(legend.position = "none") +
                       ggtitle("(c)"),
                     plot2d +
                       theme(legend.position = "none") +
                       ggtitle("(d)"),
                     plot2e +
                       theme(legend.position = "none") +
                       ggtitle("(e)"),
                     legend,
                     ncol = 2)
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figure2.jpeg"),
                plot2,
                height = 12,
                width = 10,
                bg  = "white")    


# Figure 3
## Plot
plot3 <- 
  cowplot::plot_grid(plot3a +
                       theme(legend.position = "none") +
                       ggtitle("(a)"),
                     plot3b +
                       theme(legend.position = "none") +
                       ggtitle("(b)"),
                     plot3c +
                       theme(legend.position = "none") +
                       ggtitle("(c)"),
                     plot3d +
                       theme(legend.position = "none") +
                       ggtitle("(d)"),
                     ncol = 2)
## Add legend
plot3 <- 
  cowplot::plot_grid(legend,
                     plot3,
                     nrow = 2,
                     rel_heights = c(0.1, 0.9))
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figure3.jpeg"),
                plot3,
                height = 11,
                width = 10,
                bg  = "white")  
