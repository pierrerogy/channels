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
  dplyr::mutate(larvae = as.factor(larvae),
                subsidy = as.factor(subsidy),
                shading = as.factor(shading),
                week = as.factor(week)) %>% 
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
                                   week = as.factor(week)) %>% 
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
  ## Put variables on log scale
  dplyr::mutate(across(c(din, np),
                       ~ log(.x),
                       .names = "{.col}_log"))%>% 
  ## Put variables on log scale
  dplyr::mutate(across(po4,
                       ~.x^(0.5),
                       .names = "{.col}_sqrt"))
  

# Scale data 
exp2_center <- 
  exp2 %>% 
  ## Center around 0
  dplyr::bind_cols(
    po4_scale = scale(exp2$po4_sqrt),
    nh4_scale = scale(exp2$nh4),
    no2_scale = scale(exp2$no2),
    no3_scale = scale(exp2$no3),
    din_scale = scale(exp2$din_log),
    chloro_scale = scale(exp2$chlorophyll_ugL),
    bact_scale = scale(exp2$bact),
    np_scale = scale(exp2$np_log)) %>%
  ## Remove all NAs
  dplyr::filter(!is.na(nh4_scale) & !is.na(po4_scale) &
                  !is.na(no2_scale) & !is.na(no3_scale) &
                  !is.na(bact) & !is.na(chlorophyll_ugL))

# Step 1 - Effects of nonmoz treatments on nutrients, pH and T-----------------------------
# DIN
## Fit model
din_exp2_reslight <-
  MCMCglmm::MCMCglmm(din_log ~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian",                      
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2)
## Check asssumptions
plot(din_exp2_reslight$Sol)
autocorr.diag(din_exp2_reslight$VCV)
## Test
summary(din_exp2_reslight)
## Plot
{### Data
  exp2_din_reslight_effect <- 
    as.data.frame(ggeffects::ggpredict(model = din_exp2_reslight,
                                       terms = c("shading", "subsidy"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::rename(shading = x,
                  subsidy = group)
  ### Plot
  exp2_din_reslight_plot <- 
    ggplot(data = exp2_din_reslight_effect,
           aes(x = shading, 
               y = exp(predicted)), 
           colour = subsidy) + 
    geom_point(size = 3,
               aes(colour = subsidy),
               position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = exp(conf.low), 
                      ymax = exp(conf.high),
                      colour = subsidy), 
                  width = 0.2,
                  position = position_dodge(0.5)) +
    geom_jitter(data = exp2,
                mapping = aes(x = shading,
                              y = din,
                              colour = subsidy),
                alpha = 0.3) +
    scale_y_continuous(trans = "log10") +
    ggtitle("") +
    xlab("Light exposure") +
    scale_x_discrete(labels = c("Exposed", "Shaded")) +
    ylab(expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))) +
    scale_color_manual(name = "Subsidy",
                       labels = c("Litter", "Litter + feces"), 
                       values = c("tan1", "tan4")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

# PO4
## Fit model
po4_exp2_reslight <-
  MCMCglmm::MCMCglmm(po4_sqrt ~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian",                      
                     nitt = 100000, 
                     burnin = 10000,
                     thin = 25,
                     data = exp2)
## Check asssumptions
plot(po4_exp2_reslight$Sol)
autocorr.diag(po4_exp2_reslight$VCV)
## Test
summary(po4_exp2_reslight)
## Plot
{### Data
  exp2_po4_reslight_effect <- 
    as.data.frame(ggeffects::ggpredict(po4_exp2_reslight,
                                       terms = c("shading", "subsidy"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::rename(shading = x,
                  subsidy = group)
  ### Plot
  exp2_po4_reslight_plot <- 
    ggplot(data = exp2_po4_reslight_effect,
           aes(x = shading, 
               y = predicted^2), 
           colour = subsidy) + 
    geom_point(size = 3,
               aes(colour = subsidy),
               position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = conf.low^2, ymax = conf.high^2,
                      colour = subsidy), 
                  width = 0.2,
                  position = position_dodge(0.5)) +
    geom_jitter(data = exp2,
                mapping = aes(x = shading,
                              y = po4,
                              colour = subsidy),
                alpha = 0.3) +
    ggtitle("") +
    xlab("Light exposure") +
    scale_x_discrete(labels = c("Exposed", "Shaded")) +
    ylab(expression(paste("PO"["4"]^"3-"*" concentration (", mu, "mol"*".L"^"-1"*")"))) +
    scale_color_manual(name = "Subsidy",
                       labels = c("Litter", "Litter + feces"), 
                       values = c("tan1", "tan4")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

## pH
## Fit model
pH_exp2_reslight <-
  MCMCglmm::MCMCglmm(pH~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2)
## Check asssumptions
plot(pH_exp2_reslight$Sol)
autocorr.diag(pH_exp2_reslight$VCV)
## Test
summary(pH_exp2_reslight)
## Plot
{### Data
  exp2_pH_effect_reslight <- 
    as.data.frame(ggeffects::ggpredict(pH_exp2_reslight,
                                       terms = c("shading", "subsidy"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::rename(shading = x,
                  subsidy = group)
  ### Plot
  exp2_pH_reslight_plot <- 
    ggplot(data = exp2_pH_effect_reslight,
           aes(x = shading, 
               y = predicted), 
           colour = subsidy) + 
    geom_point(size = 3,
               aes(colour = subsidy),
               position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                      colour = subsidy), 
                  width = 0.2,
                  position = position_dodge(0.5)) +
    geom_jitter(data = exp2,
                mapping = aes(x = shading,
                              y = pH,
                              colour = subsidy),
                alpha = 0.3) +
    ggtitle("") + 
    xlab("Light exposure") +
    scale_x_discrete(labels = c("Exposed", "Shaded")) +
    ylab("pH") +
    scale_color_manual(name = "Subsidy",
                       labels = c("Litter", "Litter + feces"), 
                       values = c("tan1", "tan4")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

## Temperature
## Fit model
temperature_exp2_reslight <-
  MCMCglmm::MCMCglmm(temperature_C~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 1000,
                     thin = 25,
                     data = exp2)
## Check asssumptions
plot(temperature_exp2_reslight$Sol)
autocorr.diag(temperature_exp2_reslight$VCV)
## Test
summary(temperature_exp2_reslight)
## Plot
{### Data
  exp2_temperature_effect_reslight <- 
    as.data.frame(ggeffects::ggpredict(temperature_exp2_reslight,
                                       terms = c("shading", "subsidy"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::rename(shading = x,
                  subsidy = group)
  ### Plot
  exp2_temperature_reslight_plot <- 
    ggplot(data = exp2_temperature_effect_reslight,
           aes(x = shading, 
               y = predicted), 
           colour = subsidy) + 
    geom_point(size = 3,
               aes(colour = subsidy),
               position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                      colour = subsidy), 
                  width = 0.2,
                  position = position_dodge(0.5)) +
    geom_jitter(data = exp2,
                mapping = aes(x = shading,
                              y = temperature_C,
                              colour = subsidy),
                alpha = 0.3) +
    ggtitle("") + 
    xlab("Light exposure") +
    scale_x_discrete(labels = c("Exposed", "Shaded")) +
    ylab("Temperature (°C)") +
    scale_color_manual(name = "Subsidy",
                       labels = c("Litter", "Litter + feces"), 
                       values = c("tan1", "tan4")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

## N:P ratio
hist(exp2_center$np)
## Fit model
np_exp2_reslight <-
  MCMCglmm::MCMCglmm(np_log ~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian",                      
                     nitt = 100000, 
                     burnin = 5000,
                     thin = 50,
                     data = exp2)
## Check asssumptions
plot(np_exp2_reslight$Sol)
autocorr.diag(np_exp2_reslight$VCV)
## Test
summary(np_exp2_reslight)
## Plot
{### Data
  exp2_np_reslight_effect <- 
    as.data.frame(ggeffects::ggpredict(np_exp2_reslight,
                                       terms = c("shading", "subsidy"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::rename(shading = x,
                  subsidy = group)
  ### Plot
  exp2_np_reslight_plot <- 
    ggplot(data = exp2_np_reslight_effect,
           aes(x = shading, 
               y = exp(predicted)), 
           colour = subsidy) + 
    geom_point(size = 3,
               aes(colour = subsidy),
               position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high),
                      colour = subsidy), 
                  width = 0.2,
                  position = position_dodge(0.5)) +
    geom_jitter(data = exp2,
                mapping = aes(x = shading,
                              y = np,
                              colour = subsidy),
                alpha = 0.3) +
    scale_y_continuous(trans = "log10",
                       breaks = c(0.1, 1, 5 ,15)) +
    ggtitle("") +
    xlab("Light exposure") +
    scale_x_discrete(labels = c("Exposed", "Shaded")) +
    ylab(expression(paste("N:P ratio in water"))) +
    scale_color_manual(name = "Subsidy",
                       labels = c("Litter", "Litter + feces"), 
                       values = c("tan1", "tan4")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

# Step 2.1 - Treatments on microorganisms ------------------------------------
# Chorophyll
## Fit model
chloro_exp2_treatlight <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_treatlight$Sol)
autocorr.diag(chloro_exp2_treatlight$VCV)
## Test
summary(chloro_exp2_treatlight)
## Plot
### Get coefficients
chloro_exp2_treatlight_effect <- 
  as.data.frame(ggeffects::ggpredict(chloro_exp2_treatlight,
                                     terms = c("shading", "subsidy"),
                                     type = "re",
                                     ci.level = 0.95)) %>% 
  dplyr::rename(shading = x,
                subsidy = group)
### Plot
exp2_chloro_plot_treatlight <- 
  ggplot(data = chloro_exp2_treatlight_effect,
         aes(x = shading, 
             y = predicted), 
         colour = subsidy) + 
  geom_point(size = 3,
             aes(colour = subsidy),
             position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                    colour = subsidy), 
                width = 0.2,
                position = position_dodge(0.5)) +
  geom_jitter(data = exp2,
              mapping = aes(x = shading,
                            y = chlorophyll_ugL,
                            colour = subsidy),
              alpha = 0.3) +
  ylim(0, 30) +
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Chlorophyill concentration (ug/L)") +
  scale_color_manual(name = "Subsidy",
                     labels = c("Litter", "Litter + feces"), 
                     values = c("tan1", "tan4")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


# Bacteria
## Fit model
bact_exp2_treatlight <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*subsidy,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 5000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_treatlight$Sol)
autocorr.diag(bact_exp2_treatlight$VCV)
## Test
summary(bact_exp2_treatlight)
## Plot
### Get coefficients
bact_exp2_treatlight_effect <- 
  as.data.frame(ggeffects::ggpredict(bact_exp2_treatlight,
                                     terms = c("subsidy", "shading"),
                                     type = "re",
                                     ci.level = 0.95)) %>% 
  dplyr::rename(shading = x,
                subsidy = group)
### Plot
exp2_bact_plot_treatlight <- 
  ggplot(data = bact_exp2_treatlight_effect,
         aes(x = shading, 
             y = predicted), 
         colour = subsidy) + 
  geom_point(size = 3,
             aes(colour = subsidy),
             position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                    colour = subsidy), 
                width = 0.2,
                position = position_dodge(0.5)) +
  geom_jitter(data = exp2,
              mapping = aes(x = shading,
                            y = bact,
                            colour = subsidy),
              alpha = 0.3) +
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Approx n f bacterial cells (x10e12)") +
  scale_color_manual(name = "Subsidy",
                     labels = c("Litter", "Litter + feces"), 
                     values = c("tan1", "tan4")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


# Step 2.2 - Nutrients on microorganisms ------------------------------------
# Chorophyll
## Fit model
chloro_exp2_nutlight <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*(din_scale + po4_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_nutlight$Sol)
autocorr.diag(chloro_exp2_nutlight$VCV)
## Test
summary(chloro_exp2_nutlight)

# Bacteria 
## Fit model
bact_exp2_nutlight <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*(din_scale + po4_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 10000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_nutlight$Sol)
autocorr.diag(bact_exp2_nutlight$VCV)
## Test
summary(bact_exp2_nutlight)


# Step 2.3 - Treatments and nutrients on microorganisms ------------------------------------
# Chorophyll
## Fit model
chloro_exp2_treatnutlight <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*(subsidy + din_scale + po4_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 120000, 
                     burnin = 5000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_treatnutlight$Sol)
autocorr.diag(chloro_exp2_treatnutlight$VCV)
## Test
summary(chloro_exp2_treatnutlight)


# Bacteria
## Fit model
bact_exp2_treatnutlight <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*(subsidy + din_scale + po4_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 5000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_treatnutlight$Sol)
autocorr.diag(bact_exp2_treatnutlight$VCV)
## Test
summary(bact_exp2_treatnutlight)

# Step 2.4 - Intra-microorganisms -------------------------------------------
## Fit model
chloro_exp2_bact <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*(bact_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 10000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_bact$Sol)
autocorr.diag(chloro_exp2_bact$VCV)
## Test
summary(chloro_exp2_bact)
## Plot
{### Data
  chloro_exp2_bact_effect <- 
    ggeffects::ggpredict(chloro_exp2_bact,
                         terms = c("bact_scale", "shading"),
                         type = "re",
                         ci.level = 0.95)
  ### Plot
  plot(chloro_exp2_bact_effect,
       rawdata = T) +
    ylim(0, 60) +
    ggtitle("") +
    ylab("Chlorophyll concentration (ug/L)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

## Fit model
bact_exp2_chloro <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*chloro_scale,
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_chloro$Sol)
autocorr.diag(bact_exp2_chloro$VCV)
## Test
summary(bact_exp2_chloro)
## Plot
{### Data
  bact_exp2_chloro_effect <- 
    ggeffects::ggpredict(bact_exp2_chloro,
                         terms = c("chloro_scale", "shading"),
                         type = "re",
                         ci.level = 0.95)
  ### Plot
  plot(bact_exp2_chloro_effect,
       rawdata = T) +
    ggtitle("") +
    ylim(0,40) +
    ylab("Approx n f bacterial cells (x10e12)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

# Step 2.5 - Intra-microorganisms and nutrients -------------------------------------------
# Chlorophyll
## Fit model
chloro_exp2_bactnut <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*(din_scale + po4_scale + 
                                bact_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center%>% 
                       dplyr::filter(!is.na(din_scale) & !is.na(bact)))
## Check asssumptions
plot(chloro_exp2_bactnut$Sol)
autocorr.diag(chloro_exp2_bactnut$VCV)
## Test
summary(chloro_exp2_bactnut)


# Bacteria
## Fit model
bact_exp2_chloronut <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*(din_scale + po4_scale +
                                  chloro_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 120000, 
                     burnin = 1000,
                     thin = 50,
                     data = exp2_center%>% 
                       dplyr::filter(!is.na(din_scale) & !is.na(bact)))
## Check asssumptions
plot(bact_exp2_chloronut$Sol)
autocorr.diag(bact_exp2_chloronut$VCV)
## Test
summary(bact_exp2_chloronut)

# Step 2.6 - Intra-microorganisms and treatment -------------------------------------------
# Chlorophyll
## Fit model
chloro_exp2_bacttreat <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*(subsidy + bact_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 60000, 
                     burnin = 5000,
                     thin = 25,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_bacttreat$Sol)
autocorr.diag(chloro_exp2_bacttreat$VCV)
## Test
summary(chloro_exp2_bacttreat)
## Plot
### Get coefficients
chloro_exp2_bacttreat_effect <- 
  as.data.frame(ggeffects::ggpredict(chloro_exp2_bacttreat,
                                     terms = c("bact_scale", "shading",
                                               "subsidy"),
                                     type = "re",
                                     ci.level = 0.95)) %>% 
  dplyr::rename(bact_scale = x,
                shading = group,
                subsidy = facet)
### Plot
exp2_chloro_plot_bacttreat <- 
  ggplot(data = chloro_exp2_bacttreat_effect,
         aes(x = bact_scale, 
             y = predicted), 
         colour = subsidy) + 
  geom_point(size = 3,
             aes(colour = subsidy),
             position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                    colour = subsidy), 
                width = 0.2,
                position = position_dodge(0.5)) +
  geom_jitter(data = exp2_center,
              mapping = aes(x = bact_scale,
                            y = chlorophyll_ugL,
                            colour = subsidy),
              alpha = 0.3) +
  ylim(0, 60) +
  ggtitle("") + 
  ylab("Chlorophyill concentration (ug/L)") +
  scale_color_manual(name = "Subsidy",
                     labels = c("Litter", "Litter + feces"), 
                     values = c("tan1", "tan4")) +
  facet_wrap(~shading) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# Bacteria
## Fit model
bact_exp2_chlorotreat <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*(subsidy + chloro_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 120000, 
                     burnin = 1000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_chlorotreat$Sol)
autocorr.diag(bact_exp2_chlorotreat$VCV)
## Test
summary(bact_exp2_chlorotreat)
## Plot
### Get coefficients
bact_exp2_chlorotreat_effect <- 
  as.data.frame(ggeffects::ggpredict(bact_exp2_chlorotreat,
                                     terms = c("chloro_scale", "shading",
                                               "subsidy"),
                                     type = "re",
                                     ci.level = 0.95)) %>% 
  dplyr::rename(chloro_scale = x,
                shading = group,
                subsidy = facet)
### Plot
exp2_bact_plot_chlorotreat <- 
  ggplot(data = bact_exp2_chlorotreat_effect,
         aes(x = chloro_scale, 
             y = predicted), 
         colour = subsidy) + 
  geom_point(size = 3,
             aes(colour = subsidy),
             position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                    colour = subsidy), 
                width = 0.2,
                position = position_dodge(0.5)) +
  geom_jitter(data = exp2_center,
              mapping = aes(x = chloro_scale,
                            y = bact,
                            colour = subsidy),
              alpha = 0.3) +
  ylim(0, 60) +
  ggtitle("") + 
  ylab("Chlorophyill concentration (ug/L)") +
  scale_color_manual(name = "Subsidy",
                     labels = c("Litter", "Litter + feces"), 
                     values = c("tan1", "tan4")) +
  facet_wrap(~shading) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))



# Step 2.7 - Everything on microorganisms ------------------------------------
# Chorophyll
## Fit model
chloro_exp2_all <-
  MCMCglmm::MCMCglmm(log(chlorophyll_ugL) ~
                       shading*(subsidy + din_scale + 
                                  po4_scale + bact_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 120000, 
                     burnin = 5000,
                     thin = 50,
                     data = exp2_center)
## Check asssumptions
plot(chloro_exp2_all$Sol)
autocorr.diag(chloro_exp2_all$VCV)
## Test
summary(chloro_exp2_all)

# Bacteria
## Fit model
bact_exp2_all <-
  MCMCglmm::MCMCglmm(log(bact) ~
                       shading*(subsidy + din_scale + 
                                  po4_scale + chloro_scale),
                     random = ~ week + cup_number, 
                     family = "gaussian", 
                     nitt = 240000, 
                     burnin = 10000,
                     thin = 100,
                     data = exp2_center)
## Check asssumptions
plot(bact_exp2_all$Sol)
autocorr.diag(bact_exp2_all$VCV)
## Test
summary(bact_exp2_all)





# Step 2.8 - Picking the best model ---------------------------------------
## Chlorophyll, highest to lowest
summary(chloro_exp2_nutlight)$DIC
summary(chloro_exp2_treatnutlight)$DIC
summary(chloro_exp2_treatlight)$DIC
summary(chloro_exp2_bacttreat) $DIC
summary(chloro_exp2_bact)$DIC
summary(chloro_exp2_all)$DIC
summary(chloro_exp2_bactnut)$DIC

summary(chloro_exp2_bactnut)

## Bacteria, highest to lowest
summary(bact_exp2_chloro)$DIC
summary(bact_exp2_chlorotreat)$DIC
summary(bact_exp2_treatlight)$DIC
summary(bact_exp2_nutlight)$DIC
summary(bact_exp2_chloronut)$DIC
summary(bact_exp2_treatnutlight)$DIC
summary(bact_exp2_all)$DIC

summary(bact_exp2_all)

# Step 3  - Microorganisms and nutrients on mosquitoes ----------------------------------
hist(exp2_center$n_moz)
# Model with just microorganisms
{mozmicro_model <- 
  lme4::glmer(n_moz ~ chloro_scale * bact_scale +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center)
# Assumptions
simulationOutput <- 
  DHARMa::simulateResiduals(fittedModel = mozmicro_model, 
                    n = 2000)
plot(simulationOutput,
     quantreg = F,
     rank = T)

# Tests
mozmicro_test <- 
  afex::mixed(n_moz ~ chloro_scale * bact_scale +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center,
              type = afex_options(type = "2"),
              all_fit = TRUE,
              method = "LRT")$anova_table
}

# Model with just nutrients
{moznut_model <- 
  lme4::glmer(n_moz ~ din_scale + po4_scale +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center)
# Assumptions
simulationOutput <- 
  DHARMa::simulateResiduals(fittedModel = moznut_model, 
                            n = 2000)
plot(simulationOutput,
     quantreg = F,
     rank = T)

# Tests
moznut_test <- 
  afex::mixed(n_moz ~ din_scale + po4_scale +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center,
              type = afex_options(type = "2"),
              all_fit = TRUE,
              method = "LRT")$anova_table
}

# Model with just treatment
{moztreat_model <- 
    lme4::glmer(n_moz ~ subsidy +
                  (1|week) + (1|cup_number),
                family = poisson(link = "log"),
                data = exp2_center)
  # Assumptions
  simulationOutput <- 
    DHARMa::simulateResiduals(fittedModel = moztreat_model, 
                              n = 2000)
  plot(simulationOutput,
       quantreg = F,
       rank = T)
  
  # Tests
  moztreat_test <- 
    afex::mixed(n_moz ~ subsidy +
                  (1|week) + (1|cup_number),
                family = poisson(link = "log"),
                data = exp2_center,
                type = afex_options(type = "2"),
                all_fit = TRUE,
                method = "LRT")$anova_table
}

# Model with microorganisms and light
{mozmicrolight_model <- 
  lme4::glmer(n_moz ~ shading*(chloro_scale * bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              glmerControl(optimizer = "bobyqa"),
              data = exp2_center)
# Assumptions
simulationOutput <- 
  DHARMa::simulateResiduals(fittedModel = mozmicrolight_model, 
                    n = 2000)
plot(simulationOutput,
     quantreg = F,
     rank = T)

# Tests
mozmicrolight_test <- 
  afex::mixed(n_moz ~ shading*(chloro_scale * bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center  %>% 
                dplyr::filter(!is.na(chloro_scale) &
                                !is.na(bact_scale)),
              type = afex_options(type = "2"),
              all_fit = TRUE,
              method = "LRT")$anova_table
}

# Model with nutrients and light
{moznutlight_model <- 
    lme4::glmer(n_moz ~ shading*(din_scale + po4_scale) +
                  (1|week) + (1|cup_number),
                family = poisson(link = "log"),
                data = exp2_center)
  # Assumptions
  simulationOutput <- 
    DHARMa::simulateResiduals(fittedModel = moznutlight_model, 
                              n = 2000)
  plot(simulationOutput,
       quantreg = F,
       rank = T)
  
  # Tests
  moznutlight_test <- 
    afex::mixed(n_moz ~ shading*(din_scale + po4_scale) +
                  (1|week) + (1|cup_number),
                family = poisson(link = "log"),
                data = exp2_center,
                type = afex_options(type = "2"),
                all_fit = TRUE,
                method = "LRT")$anova_table
}

# Model with treatment and light
{moztreatlight_model <- 
    lme4::glmer(n_moz ~ subsidy * shading +
                  (1|week) + (1|cup_number),
                family = poisson(link = "log"),
                data = exp2_center)
  # Assumptions
  simulationOutput <- 
    DHARMa::simulateResiduals(fittedModel = moztreatlight_model, 
                              n = 2000)
  plot(simulationOutput,
       quantreg = F,
       rank = T)
  
  # Tests
  moztreatlight_test <- 
    afex::mixed(n_moz ~ subsidy*shading +
                  (1|week) + (1|cup_number),
                family = poisson(link = "log"),
                data = exp2_center,
                type = afex_options(type = "2"),
                all_fit = TRUE,
                method = "LRT")$anova_table
}

# Model with treatment and light and nutrient
{moztreatnutlight_model <- 
    lme4::glmer(n_moz ~ (subsidy + din_scale + po4_scale) * shading +
                  (1|week) + (1|cup_number),
                family = poisson(link = "sqrt"),
                data = exp2_center)
  # Assumptions
  simulationOutput <- 
    DHARMa::simulateResiduals(fittedModel = moztreatnutlight_model, 
                              n = 2000)
  plot(simulationOutput,
       quantreg = F,
       rank = T)
  
  # Tests
  moztreatnutlight_test <- 
    afex::mixed(n_moz ~ (subsidy + din_scale + po4_scale) +
                  (1|week) + (1|cup_number),
                family = poisson(link = "sqrt"),
                data = exp2_center,
                type = afex_options(type = "2"),
                all_fit = TRUE,
                method = "LRT")$anova_table
}

# Model with microorganisms and light and nutrients
{mozmicrolightnut_model <- 
  lme4::glmer(n_moz ~ shading*(din_scale + po4_scale + 
                                 chloro_scale * bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              glmerControl(optimizer = "bobyqa"),
              data = exp2_center)
# Assumptions
simulationOutput <- 
  DHARMa::simulateResiduals(fittedModel = mozmicrolightnut_model, 
                            n = 2000)
plot(simulationOutput,
     quantreg = F,
     rank = T)

# Tests
mozmicrolightnut_test <- 
  afex::mixed(n_moz ~ shading*(din_scale + po4_scale + 
                                 chloro_scale * bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center,
              type = afex_options(type = "2"),
              all_fit = TRUE,
              method = "LRT")$anova_table
}

# Model with microorganisms and light and treatment
{mozmicrolighttreat_model <- 
  lme4::glmer(n_moz ~ shading*(subsidy + chloro_scale *
                                 bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              glmerControl(optimizer = "bobyqa"),
              data = exp2_center)
# Assumptions
simulationOutput <- 
  DHARMa::simulateResiduals(fittedModel = mozmicrolighttreat_model, 
                            n = 2000)
plot(simulationOutput,
     quantreg = F,
     rank = T)

# Tests
mozmicrolighttreat_test <- 
  afex::mixed(n_moz ~ shading*(subsidy + chloro_scale *
                                bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center,
              type = afex_options(type = "2"),
              all_fit = TRUE,
              method = "LRT")$anova_table
}

# Model with microorganisms and light and nutrients and treatment
{mozmicrolightnuttreat_model <- 
  lme4::glmer(n_moz ~ shading*(subsidy + din_scale + po4_scale + 
                                 chloro_scale * bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              glmerControl(optimizer = "bobyqa"),
              data = exp2_center)
# Assumptions
simulationOutput <- 
  DHARMa::simulateResiduals(fittedModel = mozmicrolightnuttreat_model, 
                            n = 2000)
plot(simulationOutput,
     quantreg = F,
     rank = T)

# Tests
mozmicrolightnuttreat_test <- 
  afex::mixed(n_moz ~ shading*(subsidy + din_scale + 
                                 po4_scale + chloro_scale * 
                                 bact_scale) +
                (1|week) + (1|cup_number),
              family = poisson(link = "log"),
              data = exp2_center,
              type = afex_options(type = "2"),
              all_fit = TRUE,
              method = "LRT")$anova_table
}

## Picking the best model
### Ranked from highest to lowest
AIC(moztreatnutlight_model)
AIC(mozmicrolighttreat_model)
AIC(mozmicrolight_model)
AIC(mozmicrolightnuttreat_model)
AIC(mozmicrolightnut_model)
AIC(moztreatlight_model)
AIC(mozmicro_model)
AIC(moztreat_model)
AIC(moznutlight_model)
AIC(moznut_model)




# Quick assessment of N:P ratios  --------


exp2_center %>% 
  dplyr::select(shading, subsidy, np) %>% 
  dplyr::group_by(shading, subsidy) %>% 
  dplyr::summarise_all(mean)
