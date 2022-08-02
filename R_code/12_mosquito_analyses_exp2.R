# Mosquito analyses

# Load libraries
library(tidyverse)
library(here)
library(lme4)
library(MASS)
library(ggeffects)
library(afex)
library(car)

# Load data
mozzies <- 
  read.csv(here::here("microchannels",
                      "appdata", 
                      "mosquitoes.csv")) %>% 
  ## Add average of two wings
  dplyr::mutate(wing_length = (left_wing_mm + right_wing_mm)/2)
str(mozzies)

# Time to death -----------------------------------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
hist(mozzies$time_death)
# Model
death_model <- 
  MASS::glm.nb(time_death ~ shading*subsidy,
      data = mozzies)
# Assumptions
plot(death_model)

# Tests
death_test <- 
  car::Anova(death_model)
# Plot
## Data
death_effect <- 
  ggeffects::ggeffect(death_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
death_plot <- 
  plot(death_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("") +
  ylim(0,52) +
  ## Gotta add raw data myself, jitter in ggeffects is too much
  geom_jitter(data = mozzies %>% 
                dplyr::rename(group_col = subsidy) %>% 
                dplyr::mutate(shading = ifelse(shading =="exposed",
                                               1, 2)),
              aes(x = shading,
                  y = time_death,
                  group = group_col),
              position = position_jitterdodge(dodge.width = 0.2,
                                              jitter.width = 0.1),
              size = 3,
              alpha = 0.3,
              shape = 16) +
  ylab("Age at death (days)") +
  xlab("Light exposure") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))



# Time to pupation --------------------------------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
hist(mozzies$time_pupation)
# Model
pupation_model <- 
  MASS::glm.nb(time_pupation ~ shading*subsidy,
               data = mozzies)
# Assumptions
plot(pupation_model)

# Tests
pupation_test <- 
  car::Anova(pupation_model)
# Plot
## Data
pupation_effect <- 
  ggeffects::ggeffect(pupation_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
pupation_plot <- 
  plot(pupation_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Age at pupation (days)") +
  ## Gotta add raw data myself, jitter in ggeffects is too much
  geom_jitter(data = mozzies %>% 
                dplyr::rename(group_col = subsidy) %>% 
                dplyr::mutate(shading = ifelse(shading =="exposed",
                                               1, 2)),
              aes(x = shading,
                  y = time_pupation,
                  group = group_col),
              position = position_jitterdodge(dodge.width = 0.2,
                                              jitter.width = 0.1),
              size = 3,
              alpha = 0.3,
              shape = 16) +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))

# Time to emergence -----------------------------------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
## Couldnt fit it other than that
hist(mozzies$time_emergence)
# Model
emergence_model <- 
  glm(time_emergence ~ shading*subsidy,
      data = mozzies %>% 
        dplyr::filter(!is.na(time_emergence)),
      family = poisson(link = "log"))
# Assumptions
plot(emergence_model)

# Tests
emergence_test <- 
  car::Anova(emergence_model)
# Plot
## Data
emergence_effect <- 
  ggeffects::ggeffect(emergence_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
emergence_plot <- 
  plot(emergence_effect,
       ci = T,
       rawdata = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("") +
  ylab("Age at emergence (days)") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))



# Size at death -----------------------------------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
hist(mozzies$size_mm)
# Model
sizedeath_model <- 
  lm(log(size_mm) ~ shading*subsidy*time_death,
      data = mozzies %>% 
       dplyr::filter(!is.na(size_mm)))
# Assumptions
plot(sizedeath_model)

# Tests
sizedeath_test <- 
  car::Anova(sizedeath_model)
# Plot
## Data
sizedeath_effect <- 
  ggeffects::ggeffect(sizedeath_model,
                      terms = c("time_death", "shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
pch <- 
  ifelse(mozzies$shading == "exposed",
         15, 
         19)
## Plot
sizedeath_plot <- 
  ggplot(data = mozzies,
         aes(x = time_death,
             y = jitter(size_mm, 2))) +
  geom_point(aes(
    #pch = shading,
                 colour = subsidy),
             size = 2) +
  xlab("Age at death (days)") +
  ylab("Length at death (mm)") +
  ggtitle("") + 
  ylim(1, 8) +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  # scale_shape_manual(name = "Exposition",
  #                    labels = c("Exposed", "Shaded"), 
  #                    values = c(15, 19))  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))






# Probability of dying ----------------------------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
## Make column with death yes/no
mozzies <- 
  mozzies %>% 
  dplyr::mutate(death_yn = ifelse(!is.na(time_death), 1, 0)) 
# Model
probdeath_model <- 
  glm(death_yn ~ shading*subsidy,
      family= binomial(link = "probit"),
     data = mozzies)
# Assumptions
plot(probdeath_model)
# Tests
probdeath_test <- 
  car::Anova(probdeath_model)
# Plot
## Data
probdeath_effect <- 
  ggeffects::ggeffect(probdeath_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
probdeath_plot <- 
  plot(probdeath_effect,
       ci = T,
       rawdata = F,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Percent chance of dying") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))

# Wing length --------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
## Couldnt fit it other than that
hist(log(mozzies$wing_length))
# Model
wing_model <- 
  lm(wing_length ~ sex + shading*subsidy,
        data = mozzies %>% 
          dplyr::filter(!is.na(wing_length)))
# Assumptions
plot(wing_model)

# Tests
wing_test <- 
  car::Anova(wing_model)
# Plot
## Data
wing_effect <- 
  ggeffects::ggeffect(wing_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
wing_plot <- 
  plot(wing_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Mean wing length (mm)") +
  ## Gotta add raw data myself, jitter in ggeffects is too much
  geom_jitter(data = mozzies %>% 
                dplyr::rename(group_col = subsidy) %>% 
                dplyr::mutate(shading = ifelse(shading =="exposed",
                                               1, 2)),
              aes(x = shading,
                  y = wing_length,
                  group = group_col),
              position = position_jitterdodge(dodge.width = 0.2,
                                              jitter.width = 0.1),
              size = 3,
              alpha = 0.3,
              shape = 16) +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))


# Dry mass ----------------------------------------------------------------
## Random effect not significant, boundary fit singular, so do regular GLM
## Couldnt fit it other than that
hist(log(mozzies$dry_mass_mg))
# Model
drymass_model <- 
  lm(dry_mass_mg ~ sex + shading*subsidy,
     data = mozzies %>% 
       dplyr::filter(!is.na(dry_mass_mg)))
# Assumptions
plot(drymass_model)

# Tests
drymass_test <- 
  car::Anova(drymass_model)
# Plot
## Data
drymass_effect <- 
  ggeffects::ggeffect(drymass_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozzies$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
drymass_plot <- 
  plot(drymass_effect,
       ci = T,
       rawdata = F,
       dot.size = 3,
       line.size = 1) + 
  ## Gotta add raw data myself, jitter in ggeffects is too much
  geom_jitter(data = mozzies %>% 
                dplyr::rename(group_col = subsidy) %>% 
                dplyr::mutate(shading = ifelse(shading =="exposed",
                                               1, 2)),
              aes(x = shading,
                  y = dry_mass_mg,
                  group = group_col),
              position = position_jitterdodge(dodge.width = 0.2,
                                              jitter.width = 0.1),
              size = 3,
              alpha = 0.3,
              shape = 16) +
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Adult dry mass (mg)") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))




# Pooling data per cup ---------------------------------------------------
mozpool <- 
  mozzies %>% 
  dplyr::select(cup_number, shading, subsidy,
                time_death, time_pupation, dry_mass_mg) %>% 
  ## Remove that one time where I forgot a larvae in the wrong cup
  dplyr::filter(cup_number != 46) %>% 
  ## Dummy columns to count how many things got to each point
  dplyr::mutate(death_prop = ifelse(!is.na(time_death), 1, 0),
                pup_prop = ifelse(!is.na(time_pupation), 1, 0),
                emergence_prop = ifelse(!is.na(dry_mass_mg), 1, 0),
                dry_mass_mg= ifelse(!is.na(dry_mass_mg), dry_mass_mg, 0)) %>% 
  dplyr::select(-time_death, -time_pupation) %>% 
  ## Sum everything for each cup
  dplyr::group_by(cup_number, shading, subsidy) %>% 
  dplyr::summarise_all(sum) %>% 
  dplyr::ungroup() %>% 
  ## Get proportion by dividing by 6 (initial number of larvae added)
  dplyr::mutate(death_prop = death_prop/6,
                pup_prop = pup_prop/6,
                emergence_prop = emergence_prop/6)
  
# Proportion of dying mosquitoes ------------------------------------------
# Model
prop_probdeath_model <- 
  glm(death_prop ~ shading*subsidy,
      family= binomial(link = "probit"),
      data = mozpool)
# Assumptions
plot(prop_probdeath_model)
# Tests
prop_probdeath_test <- 
  car::Anova(prop_probdeath_model)
# Plot
## Data
prop_probdeath_effect <- 
  ggeffects::ggeffect(prop_probdeath_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
prop_probdeath_plot <- 
  plot(prop_probdeath_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("") +
  ylab("Percentage of larvae that died per cup") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))

# Proportion of pupating mosquitoes -------------------------------------------------------------------------
# Model
prop_probpup_model <- 
  glm(pup_prop ~ shading*subsidy,
      family= binomial(link = "probit"),
      data = mozpool)
# Assumptions
plot(prop_probpup_model)
# Tests
prop_probpup_test <- 
  car::Anova(prop_probpup_model)
# Plot
## Data
prop_probpup_effect <- 
  ggeffects::ggeffect(prop_probpup_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
prop_probpup_plot <- 
  plot(prop_probpup_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("") +
  ylab("Percentage of larvae that pupated per cup") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))

# Proportion of emerging mosquitoes -------------------------------------------
# Model
prop_probemergence_model <- 
  glm(emergence_prop ~ shading*subsidy,
      family= binomial(link = "probit"),
      data = mozpool)
# Assumptions
plot(prop_probemergence_model)
# Tests
prop_probemergence_test <- 
  car::Anova(prop_probemergence_model)
# Plot
## Data
prop_probemergence_effect <- 
  ggeffects::ggeffect(prop_probemergence_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
prop_probemergence_plot <- 
  plot(prop_probemergence_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  xlab("") +
  ylab("Percentage of larvae that emerged per cup") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))

# Total biomass in cup -------------------------------------------------------
# Data distribution
hist(mozpool$dry_mass_mg)
# Model
tot_biomass_model <- 
  glm(dry_mass_mg+ 0.001  ~ shading*subsidy,
      family = Gamma(link = "identity"),
      data = mozpool)
# Assumptions
plot(tot_biomass_model)
# Tests
tot_biomass_test <- 
  car::Anova(tot_biomass_model)
# Plot
## Data
tot_biomass_effect <- 
  ggeffects::ggeffect(tot_biomass_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)

## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
tot_biomass_plot <- 
  plot(tot_biomass_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  ggtitle("") + 
  ## Gotta add raw data myself, jitter in ggeffects is too much
  geom_jitter(data = mozpool %>% 
               dplyr::rename(group_col = subsidy) %>% 
               dplyr::mutate(shading = ifelse(shading =="exposed",
                                              1, 2)),
             aes(x = shading,
                 y = dry_mass_mg,
                 group = group_col),
             position = position_jitterdodge(dodge.width = 0.2,
                                             jitter.width = 0.1),
             size = 3,
             alpha = 0.3,
             shape = 16) +
  xlab("Light exposure") +
  ylab("Total biomass emerged per cup (mg)") +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  scale_x_discrete(limit = c("exposed", "shaded"),
                   labels = c("Exposed", "Shaded")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))
