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

# Pooling data per cup 
mozpool <- 
  mozzies %>% 
  dplyr::select(cup_number, shading, subsidy,size_mm, wing_length,
                time_death, time_pupation, time_emergence, dry_mass_mg) %>% 
  ## Remove that one time where I forgot a larvae in the wrong cup
  dplyr::filter(cup_number != 46) %>% 
  ## Dummy columns to count how many things got to each point
  dplyr::mutate(death_prop = ifelse(!is.na(time_death), 1, 0),
                pup_prop = ifelse(!is.na(time_pupation), 1, 0),
                emergence_prop = ifelse(!is.na(dry_mass_mg), 1, 0),
                dry_mass_mg= ifelse(!is.na(dry_mass_mg), dry_mass_mg, 0)) %>% 
  ## Get proportion by dividing by 6 (initial number of larvae added)
  dplyr::mutate(death_prop = death_prop/6,
                pup_prop = pup_prop/6,
                emergence_prop = emergence_prop/6) %>% 
  ## Get overall biomass and mean time of death, pupation, emergence
  dplyr::group_by(cup_number, shading, subsidy) %>% 
  dplyr::summarise(across(.cols = dry_mass_mg, 
                          ~ sum(na.omit(.))),
                   across(.cols = c(time_death, time_pupation, 
                                    time_emergence, size_mm, wing_length),
                          ~ mean(na.omit(.)))) %>% 
  ## Replace NaN by NAs
  dplyr::mutate(across(everything(),
                       ~ifelse(is.nan(.),
                               NA, .)))

# Get proportion of emerged females to males
prop_mf <- 
  mozzies %>% 
  dplyr::filter(!is.na(sex)) %>% 
  dplyr::select(cup_number, sex) %>% 
  dplyr::group_by(cup_number, sex) %>% 
  dplyr::tally() %>%
  tidyr::pivot_wider(names_from = sex,
                     values_from = n,
                     values_fill = 0) %>% 
  dplyr::mutate(prop_f = F/(M+F),
                .keep = "none") ## removes M and F columns

# Add proportion of emerged females to males
mozpool <- 
  mozpool %>% 
  dplyr::left_join(prop_mf,
                   by = "cup_number")
  
  
# Average time to death ------------------------------------------
# Model
avg_timedeath_model <- 
  lm(log(time_death) ~ shading*subsidy,
      data = mozpool)
# Assumptions
plot(avg_timedeath_model)
# Tests
## First type 3
avg_timedeath_test <- 
  car::Anova(avg_timedeath_model,
             type = "3")
# Plot
## Data
avg_timedeath_effect <- 
  ggeffects::ggpredict(avg_timedeath_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
avg_timedeath_plot <- 
  plot(avg_timedeath_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  geom_jitter(data = mozpool%>% 
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
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Average age at death (days)") +
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


# Average size at death ------------------------------------------
# Model
avg_sizedeath_model <- 
  lm(log(size_mm) ~ shading*subsidy*time_death,
     data = mozpool)
# Assumptions
plot(avg_sizedeath_model)
# Tests
## First type 3
avg_sizedeath_test <- 
  car::Anova(avg_sizedeath_model,
             type = "3")
# Plot
## Data
avg_sizedeath_effect <- 
  ggeffects::ggpredict(avg_sizedeath_model,
                       terms = c("time_death", "subsidy"),
                       ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
avg_sizedeath_plot <-
  ggplot(data = mozpool,
         aes(x = time_death,
             y = jitter(size_mm, 2))) +
  geom_point(aes(colour = subsidy),
             size = 2) +
  xlab("Average age at death (days)") +
  ylab("Average length at death (mm)") +
  ggtitle("") + 
  ylim(1, 8) +
  scale_color_manual(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces"), 
                     values = c("tan1", "tan4")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))



# Average time to pupation ------------------------------------------
# Model
avg_timepupation_model <- 
  lm(log(time_pupation) ~ shading*subsidy,
      data = mozpool)
# Assumptions
plot(avg_timepupation_model)
# Tests
## First type 3
avg_timepupation_test <- 
  car::Anova(avg_timepupation_model,
             type = "3")
# Plot
## Data
avg_timepupation_effect <- 
  ggeffects::ggpredict(avg_timepupation_model,
                      terms = c("shading", "subsidy"),
                      ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
avg_timepupation_plot <- 
  plot(avg_timepupation_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  geom_jitter(data = mozpool%>% 
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
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Average age at pupation (days)") +
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



# Average time to emergence ------------------------------------------
# Model
avg_timeemergence_model <- 
  lm(log(time_emergence) ~ shading*subsidy,
     data = mozpool)
# Assumptions
plot(avg_timeemergence_model)
# Tests
## First type 3
avg_timeemergence_test <- 
  car::Anova(avg_timeemergence_model,
             type = "3")
# Plot
## Data
avg_timeemergence_effect <- 
  ggeffects::ggpredict(avg_timeemergence_model,
                       terms = c("shading", "subsidy"),
                       ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
avg_timeemergence_plot <- 
  plot(avg_timeemergence_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  geom_jitter(data = mozpool%>% 
                dplyr::rename(group_col = subsidy) %>% 
                dplyr::mutate(shading = ifelse(shading =="exposed",
                                               1, 2)),
              aes(x = shading,
                  y = time_emergence,
                  group = group_col),
              position = position_jitterdodge(dodge.width = 0.2,
                                              jitter.width = 0.1),
              size = 3,
              alpha = 0.3,
              shape = 16) +
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Average age at emergence (days)") +
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

# Average wing length -----------------------------------------------------
# Model
avg_winglength_model <- 
  glm(wing_length ~ prop_f + shading*subsidy ,
      family = Gamma(link = "log"),
     data = mozpool)
# Assumptions
plot(avg_winglength_model)
# Tests
avg_winglength_test <- 
  car::Anova(avg_winglength_model,
             type = "3")
# Plot
## Data
avg_winglength_effect <- 
  ggeffects::ggpredict(avg_winglength_model,
                       terms = c("shading", "subsidy"),
                       ci.level = 0.95)
## Colour
col <- 
  ifelse(mozpool$subsidy == "litter",
         "darkgreen", 
         "saddlebrown")
## Plot
avg_winglength_plot <- 
  plot(avg_winglength_effect,
       ci = T,
       dot.size = 3,
       line.size = 1) + 
  geom_jitter(data = mozpool%>% 
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
  ggtitle("") + 
  xlab("Light exposure") +
  ylab("Average length of wings (mm)") +
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

