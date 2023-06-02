# Mosquito analyses

# Load libraries
library(tidyverse)
library(here)
library(brms)
library(bayestestR)
source(here::here("R_code",
                  "functions.R"))
library(cowplot)

# Load data
mozzies <- 
  read.csv(here::here("microchannels",
                      "appdata", 
                      "mosquitoes.csv")) %>% 
  ## Add average of two wings
  dplyr::mutate(wing_length = (left_wing_mm + right_wing_mm)/2)

# Pooling data per cup 
mozclean <- 
  mozzies %>% 
  dplyr::select(cup_number, shading, subsidy,size_mm, wing_length, sex,
                time_death, time_pupation, time_emergence, dry_mass_mg) %>% 
  ## Remove that one time where I forgot a larvae in the wrong cup
  dplyr::filter(cup_number != 46) %>% 
  ## Dummy columns to count how many things got to each point
  dplyr::mutate(death = ifelse(!is.na(time_death), 0, 1),
                pup = ifelse(!is.na(time_pupation), 1, 0)) %>% 
  dplyr::rename(exposure = shading)
  
# Time to death ------------------------------------------
# Model
timedeath_model <- 
  brms::brm(log(time_death)  ~
              exposure*subsidy + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(timedeath_model)
# Check effects
bayestestR::describe_posterior(timedeath_model)
# Plot
plot4a <- 
  plot_model_nice(model = timedeath_model,
                xax = "exposure",
                yax = "time_death",
                scale = "log",
                type = "points",
                data = mozclean)

# Size at death ------------------------------------------
# Model
sizedeath_model <- 
  brms::brm(log(size_mm)  ~
              exposure*subsidy*scale(log(time_death)) + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(sizedeath_model)
# Check effects
bayestestR::describe_posterior(sizedeath_model)
# Plot
plot4b <- 
  plot_model_nice(model = sizedeath_model,
                  xax = "exposure",
                  yax = "size_mm",
                  scale = "log",
                  type = "points",
                  data = mozclean)

# Time to pupation ------------------------------------------
# Model
timepupation_model <- 
  brms::brm(log(time_pupation)  ~
              exposure*subsidy + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(timepupation_model)
# Check effects
bayestestR::describe_posterior(timepupation_model)
# Plot
plots2a <- 
  plot_model_nice(model = timepupation_model,
                  xax = "exposure",
                  yax = "time_pupation",
                  scale = "log",
                  type = "points",
                  data = mozclean)

# Average time to emergence ------------------------------------------
# Model
timeemergence_model <- 
  brms::brm(log(time_emergence)  ~
              exposure*subsidy + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(timeemergence_model)
# Check effects
bayestestR::describe_posterior(timeemergence_model)
# Plot
plots2b <- 
  plot_model_nice(model = timeemergence_model,
                  xax = "exposure",
                  yax = "time_emergence",
                  scale = "log",
                  type = "points",
                  data = mozclean)


# Probability of dying ------------------------------------------
# Model
probdeath_model <- 
  brms::brm(death  ~
              exposure*subsidy + (1|cup_number),
            iter = 2000,
            family = bernoulli(link = "logit"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(probdeath_model)
# Check effects
bayestestR::describe_posterior(probdeath_model)
# Plot
plot4c <- 
  plot_model_nice(model = probdeath_model,
                  xax = "exposure",
                  yax = "death",
                  scale = "prob",
                  type = "points",
                  data = mozclean)

# Proportion of pupating mosquitoes -------------------------------------------------------------------------
# Model
probpup_model <- 
  brms::brm(pup  ~
              exposure*subsidy + (1|cup_number),
            iter = 2000,
            family = bernoulli(link = "logit"),    
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(probpup_model)
# Check effects
bayestestR::describe_posterior(probpup_model)
# Plot
plots2c <- 
  plot_model_nice(model = probpup_model,
                  xax = "exposure",
                  yax = "pup",
                  scale = "prob",
                  type = "points",
                  data = mozclean)

# Biomass of adult -------------------------------------------------------
# Model
biomass_model <- 
  brms::brm(log(dry_mass_mg) ~
              exposure*subsidy + (1|cup_number),
            iter = 2000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(biomass_model)
# Check effects
bayestestR::describe_posterior(biomass_model)
# Plot
plot4d <- 
  plot_model_nice(model = biomass_model,
                  xax = "exposure",
                  yax = "dry_mass_mg",
                  scale = "log",
                  type = "points",
                  data = mozclean)

# Average wing length -----------------------------------------------------
# Model
wing_model <- 
  brms::brm(log(wing_length) ~
              sex + exposure*subsidy + (1|cup_number),
            iter = 4000,
            family = gaussian(link = "identity"),    
            control = list(adapt_delta = 0.97,
                           max_treedepth = 10),
            data = mozclean)
# Check assumptions
plot(wing_model)
# Check effects
bayestestR::describe_posterior(wing_model)
# Plot
plot4e <- 
  plot_model_nice(model = wing_model,
                  xax = "exposure",
                  yax = "wing_length",
                  scale = "log",
                  type = "points",
                  data = mozclean)


# Compile and save figures ------------------------------------------------
# Get legend
legend <- 
  cowplot::get_legend(plots2a)

# Figure S1
## Plot
plots2 <- 
  cowplot::plot_grid(plots2a +
                       theme(legend.position = "none") +
                       ggtitle("a"),
                     plots2b +
                       theme(legend.position = "none") +
                       ggtitle("b"),
                     plots2c +
                       theme(legend.position = "none") +
                       ggtitle("c"),
                     legend,
                     ncol = 2)
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figures2.jpeg"),
                plots2,
                height = 8,
                width = 8,
                bg  = "white")

# Figure 2
## Plot
plot4 <- 
  cowplot::plot_grid(plot4c +
                       theme(legend.position = "none") +
                       ggtitle("a"),
                     plot4a +
                       theme(legend.position = "none") +
                       ggtitle("b"),
                     plot4b +
                       theme(legend.position = "none") +
                       ggtitle("c"),
                     plot4d +
                       theme(legend.position = "none") +
                       ggtitle("d"),
                     plot4e +
                       theme(legend.position = "none") +
                       ggtitle("e"),
                     legend,
                     ncol = 2)
## Save
ggplot2::ggsave(here::here("figures",
                           "exp2_figure4.jpeg"),
                plot4,
                height = 10,
                width = 7,
                bg  = "white")    



