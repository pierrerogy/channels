# Load necessary packages
library(tidyverse)
library(ggplot2)
library(cowplot)
library(nlme)
library(DHARMa)
library(ggeffects)
library(car)
library(MASS)
library(lme4)

# Read data
taxdata <- 
  read.csv("data_exp1/bromeliad_tax_exp.csv")

# Check data
str(taxdata)
  
# Summarise data
taxdata <- 
  taxdata %>% 
  mutate(po4 = 
           ifelse(!is.na(po4_a) & !is.na(po4_b),
                  (po4_a + po4_b) / 2, 
                  ifelse(is.na(po4_a) & is.na(po4_b),
                         NA, 
                         ifelse(is.na(po4_a), 
                                po4_b, po4_a))),
         nh4 = ifelse(!is.na(nh4_a) & !is.na(nh4_b),
                      (nh4_a + nh4_b) / 2, 
                      ifelse(is.na(nh4_a) & is.na(nh4_b),
                             NA, 
                             ifelse(is.na(nh4_a), 
                                    nh4_b, nh4_a))),
         no2 = ifelse(!is.na(no2_a) & !is.na(no2_b),
                      (no2_a + no2_b) / 2, 
                      ifelse(is.na(no2_a) & is.na(no2_b),
                             NA, 
                             ifelse(is.na(no2_a), 
                                    no2_b, no2_a))),
         no3 = ifelse(!is.na(no3_a) & !is.na(no3_b),
                      (no3_a + no3_b) / 2, 
                      ifelse(is.na(no3_a) & is.na(no3_b),
                             NA, 
                             ifelse(is.na(no3_a), 
                                    no3_b, no3_a))),
         bact = ifelse(!is.na(bact_a) & !is.na(bact_b),
                       (bact_a + bact_b) / 2, 
                       ifelse(is.na(bact_a) & is.na(bact_b),
                              NA, 
                              ifelse(is.na(bact_a), 
                                     bact_b, bact_a)))) %>% 
  dplyr::select(-po4_a:-bact_b) %>% 
  mutate(subsidy_1 = ifelse(feces_1_mg == 0, 
                            "litter_only", "litter_feces"),
         subsidy_2 = ifelse(feces_2_mg == 0, 
                            "litter_only", "litter_feces"),
         subsidy_3 = ifelse(feces_3_mg == 0, 
                            "litter_only", "litter_feces")) %>% 
  unite(subsidy,subsidy_1, subsidy_2, sep = "_", remove = F)

# First - model selection and testing -------------------------------------
# Get data
first <- 
  taxdata %>% 
  filter(visit_id < 4) %>% 
  dplyr::select(-ph, -temp, - litter_3_mg, -feces_3_mg, -subsidy_3) %>% 
  na.omit() %>% 
  mutate(visit_id = as.numeric(visit_id) - 1)

# Algae
## Full model
firstmodel_chloro <- 
  lme(log(chlorophyll_uL + 0.01) ~  visit_id * (subsidy_1 + shading + vessel) + bact + nh4 + no2 + no3 + po4,
      random = ~1|bromeliad_id,
      data = first,
      method="ML")
## Model selection
firstmodel_chloro <- stepAIC(firstmodel_chloro, 
                             direction="backward")
## Assumptions
plot(firstmodel_chloro)
## Summary and tests
summary(firstmodel_chloro)
Anova(firstmodel_chloro)


# Bacteria
## Full model
firstmodel_bact <- 
  lme(log(bact + 1) ~ visit_id * (subsidy_1 + shading + vessel) + chlorophyll_uL + nh4 + no2 + no3 + po4,
      random = ~1|bromeliad_id,
      data = first,
      method="ML")
## Model selection
firstmodel_bact <- stepAIC(firstmodel_bact, 
                             direction="backward")
## Assumptions
plot(firstmodel_bact)
## Summary and tests
summary(firstmodel_bact)
Anova(firstmodel_bact)

# Ammonium
## Full model
firstmodel_nh4 <- 
  lme(log(nh4 + 0.01) ~  visit_id * (subsidy_1 + shading + vessel) + chlorophyll_uL + bact + no2 + no3 + po4,
      random = ~1|bromeliad_id,
      data = first,
      method="ML")
## Model selection
firstmodel_nh4 <- stepAIC(firstmodel_nh4, 
                           direction="backward")
## Assumptions
plot(firstmodel_nh4)
## Summary and tests
summary(firstmodel_nh4)
Anova(firstmodel_nh4)

# Phosphate
## Full model
firstmodel_po4 <- 
  lme(log(po4 + 0.01) ~  visit_id * (subsidy_1 + shading + vessel) + chlorophyll_uL + bact + no2 + no3 + nh4,
      random = ~1|bromeliad_id,
      data = first,
      method="ML")
## Model selection
firstmodel_po4 <- stepAIC(firstmodel_po4, 
                          direction="backward")
## Assumptions
plot(firstmodel_po4)
## Summary and tests
summary(firstmodel_po4)
Anova(firstmodel_po4)

# Nitrite
## Full model
firstmodel_no2 <- 
  lme(log(no2 + 0.01) ~  visit_id * (subsidy_1 + shading + vessel) + chlorophyll_uL + bact + po4 + no3 + nh4,
      random = ~1|bromeliad_id,
      data = first,
      method="ML")
## Model selection
firstmodel_no2 <- stepAIC(firstmodel_no2, 
                          direction="backward")
## Assumptions
plot(firstmodel_no2)
## Summary and tests
summary(firstmodel_no2)
Anova(firstmodel_no2)

# Nitrate
## Full model
firstmodel_no3 <- 
  lme(no3  ~  visit_id * (subsidy_1 + shading + vessel) + chlorophyll_uL + bact + po4 + no2 + nh4,
      random = ~1|bromeliad_id,
      data = first,
      method="ML")
## Model selection
firstmodel_no3 <- stepAIC(firstmodel_no3, 
                          direction="backward")
## Assumptions
plot(firstmodel_no3)
## Summary and tests
summary(firstmodel_no3)
Anova(firstmodel_no3)





# Second - model selection and testing ------------------------------------
# Get data
second <- 
  taxdata %>% 
  filter(visit_id > 3) %>% 
  dplyr::select(-ph, -temp, - litter_3_mg, -feces_3_mg, -subsidy_3) %>% 
  na.omit() %>% 
  mutate(visit_id = as.numeric(visit_id) - 1)

# Algae
## Full model
secondmodel_chloro <- 
  lme(log(chlorophyll_uL + 0.01) ~  visit_id * (subsidy + shading) + bact + nh4 + no2 + no3 + po4,
      random = ~1|bromeliad_id,
      data = second,
      method="ML")
## Model selection
secondmodel_chloro <- stepAIC(secondmodel_chloro, 
                             direction="backward")
## Assumptions
plot(secondmodel_chloro)
## Summary and tests
summary(secondmodel_chloro)
Anova(secondmodel_chloro)

# Bacteria
## Full model
secondmodel_bact <- 
  lme(log(bact + 1) ~  visit_id * (subsidy + shading) + chlorophyll_uL + nh4 + no2 + no3 + po4,
      random = ~1|bromeliad_id,
      data = second,
      method="ML")
## Model selection
secondmodel_bact <- stepAIC(secondmodel_bact, 
                              direction="backward")
## Assumptions
plot(secondmodel_bact)
## Summary and tests
summary(secondmodel_bact)
Anova(secondmodel_bact)

# Ammonium
## Full model
secondmodel_nh4 <- 
  lme(log(nh4 + 0.01) ~  visit_id * (subsidy + shading) + chlorophyll_uL + bact + no2 + no3 + po4,
      random = ~1|bromeliad_id,
      data = second,
      method="ML")
## Model selection
secondmodel_nh4 <- stepAIC(secondmodel_nh4, 
                            direction="backward")
## Assumptions
plot(secondmodel_nh4)
## Summary and tests
summary(secondmodel_nh4)
Anova(secondmodel_nh4)

# Phosphate
## Full model
secondmodel_po4 <- 
  lme(log(po4 + 0.01) ~  visit_id * (subsidy + shading) + chlorophyll_uL + bact + no2 + no3 + nh4,
      random = ~1|bromeliad_id,
      data = second,
      method="ML")
## Model selection
secondmodel_po4 <- stepAIC(secondmodel_po4, 
                           direction="backward")
## Assumptions
plot(secondmodel_po4)
## Summary and tests
summary(secondmodel_po4)
Anova(secondmodel_po4)

# Nitrite
## Full model
secondmodel_no2 <- 
  lme(log(no2 + 0.01) ~  visit_id * (subsidy + shading) + chlorophyll_uL + bact + po4 + no3 + nh4,
      random = ~1|bromeliad_id,
      data = second,
      method="ML")
## Model selection
secondmodel_no2 <- stepAIC(secondmodel_no2, 
                           direction="backward")
## Assumptions
plot(secondmodel_no2)
## Summary and tests
summary(secondmodel_no2)
Anova(secondmodel_no2)

# Nitrate
## Full model
secondmodel_no3 <- 
  lme(log(no3 + 0.01) ~  visit_id * (subsidy + shading) + chlorophyll_uL + bact + po4 + no2 + nh4,
      random = ~1|bromeliad_id,
      data = second,
      method="ML")
## Model selection
secondmodel_no3 <- stepAIC(secondmodel_no3, 
                           direction="backward")
## Assumptions
plot(secondmodel_no3)
## Summary and tests
summary(secondmodel_no3)
Anova(secondmodel_no3)


# Plots - Algae -------------------------------------------------------------------
# Visit_id and shading
## First
###Data
chloroeffect_first_visitid_shading <- 
  ggeffect(firstmodel_chloro,
           terms = c("visit_id", "shading"),
           type = "re",
           ci.level = 0.95)
chloroeffect_first_visitid_shading$conf.low <- 
  exp(chloroeffect_first_visitid_shading$conf.low) 
chloroeffect_first_visitid_shading$conf.high <- 
  exp(chloroeffect_first_visitid_shading$conf.high)
chloroeffect_first_visitid_shading$predicted <- 
  exp(chloroeffect_first_visitid_shading$predicted)
### Plot
chloroplot_first_visitid_shading <- 
  plot(chloroeffect_first_visitid_shading,
       ci = T) + 
  geom_jitter(data = first,
             mapping = aes(x = as.numeric(visit_id), 
                           y = chlorophyll_uL,
                           col = shading,
                           fill = shading)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  scale_x_continuous(breaks = c(1, 2)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Chlorophyll (ug.L)") +
  scale_color_manual(name ="Shading",
                     labels = c("Light", "Shade"), 
                     values = c("goldenrod", "grey50")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
## Second
### Data
chloroeffect_second_visitid_shading <- 
  ggeffect(secondmodel_chloro,
           terms = c("visit_id", "shading"),
           type = "re",
           ci.level = 0.95)
chloroeffect_second_visitid_shading$conf.low <- 
  exp(chloroeffect_second_visitid_shading$conf.low) 
chloroeffect_second_visitid_shading$conf.high <- 
  exp(chloroeffect_second_visitid_shading$conf.high)
chloroeffect_second_visitid_shading$predicted <- 
  exp(chloroeffect_second_visitid_shading$predicted)
### Plot
chloroplot_second_visitid_shading <- 
  plot(chloroeffect_second_visitid_shading,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = as.numeric(visit_id), 
                            y = chlorophyll_uL,
                            col = shading,
                            fill = shading)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30, 100)) +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  ggtitle("Second") + 
  xlab("Week") +
  ylab("Chlorophyll (ug.L)") +
  scale_color_manual(name ="Shading",
                     labels = c("Light", "Shade"), 
                     values = c("goldenrod", "grey50")) +
  theme(legend.position = c(0.9,1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Vessel and subsidy
## First
### Data
chloroeffect_first_vessel_subsidy <- 
  ggeffect(firstmodel_chloro,
           terms = c("vessel", "subsidy_1"),
           type = "re",
           ci.level = 0.95)
chloroeffect_first_vessel_subsidy$conf.low <- 
  exp(chloroeffect_first_vessel_subsidy$conf.low) 
chloroeffect_first_vessel_subsidy$conf.high <- 
  exp(chloroeffect_first_vessel_subsidy$conf.high)
chloroeffect_first_vessel_subsidy$predicted <- 
  exp(chloroeffect_first_vessel_subsidy$predicted)
### Plot
chloroplot_first_vessel_subsidy <- 
  plot(chloroeffect_first_vessel_subsidy,
       ci = T) + 
  scale_x_discrete(labels = c("Bromeliad",
                              "Wax")) +
  geom_jitter(data = first,
              mapping = aes(x = vessel, 
                            y = chlorophyll_uL,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  
  ggtitle("First") + 
  xlab("Vessel") +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  ylab("Chlorophyll (ug.L)") +
  theme(legend.position = c(0.9,1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
## Second
### Data
chloroeffect_second_vessel_subsidy <- 
  ggeffect(secondmodel_chloro,
           terms = c("subsidy"),
           type = "re",
           ci.level = 0.95)
chloroeffect_second_vessel_subsidy$conf.low <- 
  exp(chloroeffect_second_vessel_subsidy$conf.low) 
chloroeffect_second_vessel_subsidy$conf.high <- 
  exp(chloroeffect_second_vessel_subsidy$conf.high)
chloroeffect_second_vessel_subsidy$predicted <- 
  exp(chloroeffect_second_vessel_subsidy$predicted)
### Plot
chloroplot_second_vessel_subsidy <- 
  plot(chloroeffect_second_vessel_subsidy,
       ci = T) + 
  scale_x_discrete(labels = c("Litter + Feces/Litter + Feces ",
                              "Litter + Feces/Litter",
                              "Litter/Litter + Feces",
                              "Litter/Litter")) +
  geom_jitter(data = second,
              mapping = aes(x = subsidy, 
                            y = chlorophyll_uL,
                            col = subsidy,
                            fill = subsidy)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5, 30, 90)) +
  scale_colour_manual(values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) +
  ggtitle("Second") + 
  xlab("Subsidy 1/Subsidy 2") +
  ylab("Chlorophyll (ug.L)") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Bacteria
## First
### Data
chloroeffect_first_bacteria <- 
  ggeffect(firstmodel_chloro,
           terms = "bact",
           type = "re",
           ci.level = 0.95)
chloroeffect_first_bacteria$conf.low <- 
  exp(chloroeffect_first_bacteria$conf.low) 
chloroeffect_first_bacteria$conf.high <- 
  exp(chloroeffect_first_bacteria$conf.high)
chloroeffect_first_bacteria$predicted <- 
  exp(chloroeffect_first_bacteria$predicted)
### Plot
chloroplot_first_bacteria <- 
  plot(chloroeffect_first_bacteria,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = bact, 
                            y = chlorophyll_uL,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("First") + 
  xlab("Estimated bacteria cell density (/10e11)") +
  ylab("Chlorophyll (ug.L)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
## Second
### Data
chloroeffect_second_bacteria <- 
  ggeffect(secondmodel_chloro,
           terms = c("bact"),
           type = "re",
           ci.level = 0.95)
chloroeffect_second_bacteria$conf.low <- 
  exp(chloroeffect_second_bacteria$conf.low) 
chloroeffect_second_bacteria$conf.high <- 
  exp(chloroeffect_second_bacteria$conf.high)
chloroeffect_second_bacteria$predicted <- 
  exp(chloroeffect_second_bacteria$predicted)
### Plot
chloroplot_second_bacteria <- 
  plot(chloroeffect_second_bacteria,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = bact, 
                            y = chlorophyll_uL,
                            col = shading,
                            fill = shading)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("Second") + 
  xlab("Estimated bacteria cell density (/10e11)") +
  ylab("Chlorophyll (ug.L)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Nitrite
## First
### Data
chloroeffect_first_nitrite <- 
  ggeffect(firstmodel_chloro,
           terms = "no2",
           type = "re",
           ci.level = 0.95)
chloroeffect_first_nitrite$conf.low <- 
  exp(chloroeffect_first_nitrite$conf.low) 
chloroeffect_first_nitrite$conf.high <- 
  exp(chloroeffect_first_nitrite$conf.high)
chloroeffect_first_nitrite$predicted <- 
  exp(chloroeffect_first_nitrite$predicted)
### Plot
chloroplot_first_nitrite <- 
  plot(chloroeffect_first_nitrite,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = no2, 
                            y = chlorophyll_uL,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("First") + 
  xlab("Nitrite concentration (umol.L)") +
  ylab("Chlorophyll (ug.L)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
## Second
### Data
chloroeffect_second_nitrite <- 
  ggeffect(secondmodel_chloro,
           terms = c("no2"),
           type = "re",
           ci.level = 0.95)
chloroeffect_second_nitrite$conf.low <- 
  exp(chloroeffect_second_nitrite$conf.low) 
chloroeffect_second_nitrite$conf.high <- 
  exp(chloroeffect_second_nitrite$conf.high)
chloroeffect_second_nitrite$predicted <- 
  exp(chloroeffect_second_nitrite$predicted)
### Plot
chloroplot_second_nitrite <- 
  plot(chloroeffect_second_nitrite,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = no2, 
                            y = chlorophyll_uL,
                            col = shading,
                            fill = shading)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("Second") + 
  xlab("Nitrite concentration (umol.L)") +
  ylab("Chlorophyll (ug.L)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))





# Plots - Bacteria --------------------------------------------------------
# Visit_id and vessel
## First
###Data
bacteffect_first_visitid_vessel <- 
  ggeffect(firstmodel_bact,
           terms = c("visit_id", "vessel"),
           type = "re",
           ci.level = 0.95)
bacteffect_first_visitid_vessel$conf.low <- 
  exp(bacteffect_first_visitid_vessel$conf.low) 
bacteffect_first_visitid_vessel$conf.high <- 
  exp(bacteffect_first_visitid_vessel$conf.high)
bacteffect_first_visitid_vessel$predicted <- 
  exp(bacteffect_first_visitid_vessel$predicted)
### Plot
bactplot_first_visitid_vessel <- 
  plot(bacteffect_first_visitid_vessel,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = as.numeric(visit_id), 
                            y = bact,
                            col = vessel,
                            fill = vessel)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5, 10, 600)) +
  scale_x_continuous(breaks = c(1, 2)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_color_manual(name ="Vessel",
                     labels = c("Bromeliad", "Wax"), 
                     values = c("darkgreen", "navajowhite2")) +
  theme(legend.position = c(0.9, 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Subsidy
## First
### Data
bacteffect_first_subsidy <- 
  ggeffect(firstmodel_bact,
           terms = c("subsidy_1"),
           type = "re",
           ci.level = 0.95)
bacteffect_first_subsidy$conf.low <- 
  exp(bacteffect_first_subsidy$conf.low) 
bacteffect_first_subsidy$conf.high <- 
  exp(bacteffect_first_subsidy$conf.high)
bacteffect_first_subsidy$predicted <- 
  exp(bacteffect_first_subsidy$predicted)
### Plot
bactplot_first_subsidy <- 
  plot(bacteffect_first_subsidy,
       ci = T) + 
  scale_x_discrete(labels = c("Litter + Feces",
                              "Litter")) +
  geom_jitter(data = first,
              mapping = aes(x = subsidy_1, 
                            y = bact,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(100, 200, 600)) +
  ggtitle("First") + 
  xlab("Subsidy 1") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_color_manual(values = c("burlywood", "cadetblue")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
## Second
### Data
bacteffect_second_subsidy <- 
  ggeffect(secondmodel_bact,
           terms = c("subsidy"),
           type = "re",
           ci.level = 0.95)
bacteffect_second_subsidy$conf.low <- 
  exp(bacteffect_second_subsidy$conf.low) 
bacteffect_second_subsidy$conf.high <- 
  exp(bacteffect_second_subsidy$conf.high)
bacteffect_second_subsidy$predicted <- 
  exp(bacteffect_second_subsidy$predicted)
### Plot
bactplot_second_subsidy <- 
  plot(bacteffect_second_subsidy,
       ci = T) + 
  scale_x_discrete(labels = c("Litter + Feces/Litter + Feces ",
                              "Litter + Feces/Litter",
                              "Litter/Litter + Feces",
                              "Litter/Litter")) +
  geom_jitter(data = second,
              mapping = aes(x = subsidy, 
                            y = bact,
                            col = subsidy,
                            fill = subsidy)) +
  scale_y_continuous(trans = "log",
                     breaks = c(100, 200, 600)) +
  ggtitle("Second") + 
  xlab("Subsidy 1/Subsidy 2") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_color_manual(values = c("burlywood", "grey50",
                                "cadetblue2", "cadetblue")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Ammonium
## First
### Data
bacteffect_first_ammonium <- 
  ggeffect(firstmodel_bact,
           terms = "nh4",
           type = "re",
           ci.level = 0.95)
bacteffect_first_ammonium$conf.low <- 
  exp(bacteffect_first_ammonium$conf.low) 
bacteffect_first_ammonium$conf.high <- 
  exp(bacteffect_first_ammonium$conf.high)
bacteffect_first_ammonium$predicted <- 
  exp(bacteffect_first_ammonium$predicted)
### Plot
bactplot_first_ammonium <- 
  plot(bacteffect_first_ammonium,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = nh4, 
                            y = bact,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("") + 
  xlab("Ammonium concentration (umol.L)") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Nitrite
## Second
### Data
bacteffect_second_nitrite <- 
  ggeffect(secondmodel_bact,
           terms = "no2",
           type = "re",
           ci.level = 0.95)
bacteffect_second_nitrite$conf.low <- 
  exp(bacteffect_second_nitrite$conf.low) 
bacteffect_second_nitrite$conf.high <- 
  exp(bacteffect_second_nitrite$conf.high)
bacteffect_second_nitrite$predicted <- 
  exp(bacteffect_second_nitrite$predicted)
### Plot
bactplot_second_nitrite <- 
  plot(bacteffect_second_nitrite,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = no2, 
                            y = bact,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("") + 
  xlab("Nitrite concentration (umol.L)") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Nitrate
## First
### Data
bacteffect_first_nitrate <- 
  ggeffect(firstmodel_bact,
           terms = "no3",
           type = "re",
           ci.level = 0.95)
bacteffect_first_nitrate$conf.low <- 
  exp(bacteffect_first_nitrate$conf.low) 
bacteffect_first_nitrate$conf.high <- 
  exp(bacteffect_first_nitrate$conf.high)
bacteffect_first_nitrate$predicted <- 
  exp(bacteffect_first_nitrate$predicted)
### Plot
bactplot_first_nitrate <- 
  plot(bacteffect_first_nitrate,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = no3, 
                            y = bact,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("") + 
  xlab("Nitrate concentration (umol.L)") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
## Second
### Data
bacteffect_second_nitrate <- 
  ggeffect(secondmodel_bact,
           terms = c("no3"),
           type = "re",
           ci.level = 0.95)
bacteffect_second_nitrate$conf.low <- 
  exp(bacteffect_second_nitrate$conf.low) 
bacteffect_second_nitrate$conf.high <- 
  exp(bacteffect_second_nitrate$conf.high)
bacteffect_second_nitrate$predicted <- 
  exp(bacteffect_second_nitrate$predicted)
### Plot
bactplot_second_nitrate <- 
  plot(bacteffect_second_nitrate,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = no3, 
                            y = bact,
                            col = shading,
                            fill = shading)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,30,100)) +
  ggtitle("") + 
  xlab("Nitrate concentration (umol.L)") +
  ylab("Estimated bacteria cell density (/10e11)") +
  scale_colour_manual(values = c("grey50", "grey50")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Plots - Ammonium and treatment effects ---------------------------------
# Visit id and vessel
## First
### Data
ammoniumeffect_first_visitid_vessel <- 
  ggeffect(firstmodel_nh4,
           terms = c("visit_id", "vessel"),
           type = "re",
           ci.level = 0.95)
ammoniumeffect_first_visitid_vessel$conf.low <- 
  exp(ammoniumeffect_first_visitid_vessel$conf.low) 
ammoniumeffect_first_visitid_vessel$conf.high <- 
  exp(ammoniumeffect_first_visitid_vessel$conf.high)
ammoniumeffect_first_visitid_vessel$predicted <- 
  exp(ammoniumeffect_first_visitid_vessel$predicted)
### Plot
ammoniumplot_first_visitid_vessel <- 
  plot(ammoniumeffect_first_visitid_vessel,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = as.numeric(visit_id), 
                            y = nh4,
                            col = vessel,
                            fill = vessel)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Ammonium concentratin (umol.L)") +
  scale_x_continuous(breaks = c(1, 2)) +
  scale_color_manual(name ="Vessel",
                     labels = c("Bromeliad", "Wax"), 
                     values = c("darkgreen", "navajowhite2")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


# Visit id and subsidy
### First
### Data
ammoniumeffect_first_visitid_subsidy <- 
  ggeffect(firstmodel_nh4,
           terms = c("visit_id", "subsidy_1"),
           type = "re",
           ci.level = 0.95)
ammoniumeffect_first_visitid_subsidy$conf.low <- 
  exp(ammoniumeffect_first_visitid_subsidy$conf.low) 
ammoniumeffect_first_visitid_subsidy$conf.high <- 
  exp(ammoniumeffect_first_visitid_subsidy$conf.high)
ammoniumeffect_first_visitid_subsidy$predicted <- 
  exp(ammoniumeffect_first_visitid_subsidy$predicted)
### Plot
ammoniumplot_first_visitid_subsidy <- 
  plot(ammoniumeffect_first_visitid_subsidy,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = as.numeric(visit_id), 
                            y = nh4,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Ammonium concentratin (umol.L)") +
  scale_color_manual(name ="Subsidy 1",
                     labels = c("Litter + Feces", "Litter"), 
                     values = c("burlywood", "cadetblue")) +
  scale_x_continuous(breaks = c(1, 2)) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
### Second
### Data
ammoniumeffect_second_visitid_subsidy <- 
  ggeffect(secondmodel_nh4,
           terms = c("visit_id"),
           type = "re",
           ci.level = 0.95)
ammoniumeffect_second_visitid_subsidy$conf.low <- 
  exp(ammoniumeffect_second_visitid_subsidy$conf.low) 
ammoniumeffect_second_visitid_subsidy$conf.high <- 
  exp(ammoniumeffect_second_visitid_subsidy$conf.high)
ammoniumeffect_second_visitid_subsidy$predicted <- 
  exp(ammoniumeffect_second_visitid_subsidy$predicted)
### Plot
ammoniumplot_second_visitid_subsidy <- 
  plot(ammoniumeffect_second_visitid_subsidy,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = as.numeric(visit_id), 
                            y = nh4,
                            col = subsidy,
                            fill = subsidy)) +
  ggtitle("Second") + 
  xlab("Week") +
  ylab("Ammonium concentratin (umol.L)") +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))




# Plots - Phosphate and treatment effects ---------------------------------
# Visit id and vessel
## First
### Data
phosphateeffect_first_visitid_vessel <- 
  ggeffect(firstmodel_po4,
           terms = c("visit_id", "vessel"),
           type = "re",
           ci.level = 0.95)
phosphateeffect_first_visitid_vessel$conf.low <- 
  exp(phosphateeffect_first_visitid_vessel$conf.low) 
phosphateeffect_first_visitid_vessel$conf.high <- 
  exp(phosphateeffect_first_visitid_vessel$conf.high)
phosphateeffect_first_visitid_vessel$predicted <- 
  exp(phosphateeffect_first_visitid_vessel$predicted)
### Plot
phosphateplot_first_visitid_vessel <- 
  plot(phosphateeffect_first_visitid_vessel,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = as.numeric(visit_id), 
                            y = po4,
                            col = vessel,
                            fill = vessel)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Phosphate concentratin (umol.L)") +
  scale_x_continuous(breaks = c(1, 2)) +
  scale_color_manual(name ="Vessel",
                     labels = c("Bromeliad", "Wax"), 
                     values = c("darkgreen", "navajowhite2")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


# Visit id and subsidy
### First
### Data
phosphateeffect_first_visitid_subsidy <- 
  ggeffect(firstmodel_po4,
           terms = c("visit_id", "subsidy_1"),
           type = "re",
           ci.level = 0.95)
phosphateeffect_first_visitid_subsidy$conf.low <- 
  exp(phosphateeffect_first_visitid_subsidy$conf.low) 
phosphateeffect_first_visitid_subsidy$conf.high <- 
  exp(phosphateeffect_first_visitid_subsidy$conf.high)
phosphateeffect_first_visitid_subsidy$predicted <- 
  exp(phosphateeffect_first_visitid_subsidy$predicted)
### Plot
phosphateplot_first_visitid_subsidy <- 
  plot(phosphateeffect_first_visitid_subsidy,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = as.numeric(visit_id), 
                            y = po4,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Phosphate concentratin (umol.L)") +
  scale_x_continuous(breaks = c(1, 2)) +
  scale_color_manual(name ="Subsidy 1",
                     labels = c("Litter + Feces", "Litter"), 
                     values = c("burlywood", "cadetblue")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
### Second
### Data
phosphateeffect_second_visitid_subsidy <- 
  ggeffect(secondmodel_po4,
           terms = c("visit_id", "subsidy"),
           type = "re",
           ci.level = 0.95)
phosphateeffect_second_visitid_subsidy$conf.low <- 
  exp(phosphateeffect_second_visitid_subsidy$conf.low) 
phosphateeffect_second_visitid_subsidy$conf.high <- 
  exp(phosphateeffect_second_visitid_subsidy$conf.high)
phosphateeffect_second_visitid_subsidy$predicted <- 
  exp(phosphateeffect_second_visitid_subsidy$predicted)
### Plot
phosphateplot_second_visitid_subsidy <- 
  plot(phosphateeffect_second_visitid_subsidy,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = as.numeric(visit_id), 
                            y = po4,
                            col = subsidy,
                            fill = subsidy)) +
  ggtitle("Second") + 
  xlab("Week") +
  ylab("Phosphate concentratin (umol.L)") +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))





# Plots - Nitrite and treatment effects ---------------------------------
# Visit id and subsidy
### First
### Data
nitriteeffect_first_visitid_subsidy <- 
  ggeffect(firstmodel_no2,
           terms = c("visit_id", "subsidy_1"),
           type = "re",
           ci.level = 0.95)
nitriteeffect_first_visitid_subsidy$conf.low <- 
  exp(nitriteeffect_first_visitid_subsidy$conf.low) 
nitriteeffect_first_visitid_subsidy$conf.high <- 
  exp(nitriteeffect_first_visitid_subsidy$conf.high)
nitriteeffect_first_visitid_subsidy$predicted <- 
  exp(nitriteeffect_first_visitid_subsidy$predicted)
### Plot
nitriteplot_first_visitid_subsidy <- 
  plot(nitriteeffect_first_visitid_subsidy,
       ci = T) + 
  geom_jitter(data = first,
              mapping = aes(x = as.numeric(visit_id), 
                            y = no2,
                            col = subsidy_1,
                            fill = subsidy_1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,10,15)) +
  scale_x_continuous(breaks = c(1, 2)) +
  ggtitle("First") + 
  xlab("Week") +
  ylab("Nitrite concentratin (umol.L)") +
  scale_color_manual(name ="Subsidy 1",
                     labels = c("Litter + Feces", "Litter"), 
                     values = c("burlywood", "cadetblue")) +
  theme(legend.position = c(0.9, 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
### Second
### Data
nitriteeffect_second_visitid_subsidy <- 
  ggeffect(secondmodel_no2,
           terms = c("visit_id", "subsidy"),
           type = "re",
           ci.level = 0.95)
nitriteeffect_second_visitid_subsidy$conf.low <- 
  exp(nitriteeffect_second_visitid_subsidy$conf.low) 
nitriteeffect_second_visitid_subsidy$conf.high <- 
  exp(nitriteeffect_second_visitid_subsidy$conf.high)
nitriteeffect_second_visitid_subsidy$predicted <- 
  exp(nitriteeffect_second_visitid_subsidy$predicted)
### Plot
nitriteplot_second_visitid_subsidy <- 
  plot(nitriteeffect_second_visitid_subsidy,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = as.numeric(visit_id), 
                            y = no2,
                            col = subsidy,
                            fill = subsidy)) +
  ggtitle("Second") + 
  xlab("Week") +
  ylab("Nitrite concentratin (umol.L)") +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# Visit id and shading
## Second
### Data
nitriteeffect_second_visitid_shading <- 
  ggeffect(secondmodel_no2,
           terms = c("visit_id", "shading"),
           type = "re",
           ci.level = 0.95)
nitriteeffect_second_visitid_shading$conf.low <- 
  exp(nitriteeffect_second_visitid_shading$conf.low) 
nitriteeffect_second_visitid_shading$conf.high <- 
  exp(nitriteeffect_second_visitid_shading$conf.high)
nitriteeffect_second_visitid_shading$predicted <- 
  exp(nitriteeffect_second_visitid_shading$predicted)
### Plot
nitriteplot_second_visitid_shading <- 
  plot(nitriteeffect_second_visitid_shading,
       ci = T) + 
  geom_jitter(data = second,
              mapping = aes(x = as.numeric(visit_id), 
                            y = no2,
                            col = shading,
                            fill = shading)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5,10,15)) +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  ggtitle("Second") + 
  xlab("Week") +
  ylab("Nitrite concentratin (umol.L)") +
  scale_color_manual(name ="shading",
                     labels = c("Bromeliad", "Wax"), 
                     values = c("darkgreen", "navajowhite2")) +
  theme(legend.position = c(0.9, 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))





# Plots - Associations between nutrients ----------------------------------
# Ammonium
plot_second_ammonium_nitrite <- 
  ggplot(second, 
         aes(x = nh4, 
             y = no2,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  xlab("Ammonium concentration (umol.L)") +
  ylab("Nitrite concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw() 

# Phosphate
plot_first_phosphate_ammonium <- 
  ggplot(first, 
         aes(x = po4, 
             y = nh4,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  ylab("Ammonium concentration (umol.L)") +
  xlab("Phosphate concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw() 

plot_second_phosphate_nitrite <- 
  ggplot(second, 
         aes(x = po4, 
             y = no2,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  ylab("Nitrite concentration (umol.L)") +
  xlab("Phosphate concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw() 

plot_second_phosphate_nitrate <- 
  ggplot(second, 
         aes(x = po4, 
             y = no3,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  ylab("Nitrate concentration (umol.L)") +
  xlab("Phosphate concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw() 

# Nitrite
plot_first_nitrite_ammonium <- 
  ggplot(first, 
         aes(x = no2, 
             y = nh4,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  ylab("Ammoniun concentration (umol.L)") +
  xlab("Nitrite concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw()

plot_second_nitrite_ammonium <- 
  ggplot(second, 
         aes(x = no2, 
             y = nh4,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  ylab("Ammoniun concentration (umol.L)") +
  xlab("Nitrite concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw()

plot_second_nitrite_phosphate <- 
  ggplot(second, 
         aes(x = no2, 
             y = po4,
             shape = vessel,
             colour = subsidy,
             group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) + 
  scale_colour_manual(name = "Subsidy 1/ Subsidy 2",
                      labels = c("Litter + Feces/Litter + Feces",  "Litter + Feces/Litter",
                                 "Litter/Litter + Feces", "Litter/Litter"),
                      values = c("burlywood", "grey50",
                                 "cadetblue2", "cadetblue")) + 
  ylab("Phosphate concentration (umol.L)") +
  xlab("Nitrite concentration (umol.L)") +
  scale_x_discrete() +
  theme_bw()



# Plot time series --------------------------------------------------------
# Only filter corresponding data, convert to numeric and remove unused columns
time_series <- 
  taxdata %>%
  filter(grepl("H", visit_id)) %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H",
            replacement = "") %>%
  mutate(visit_id = as.numeric(visit_id)) %>% 
  dplyr::select(-litter_1_mg:-feces_2_mg)

# Plot pH
timeseries_ph <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = ph,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("pH") +
  theme_bw() 

# Plot temperature
timeseries_temp <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = temp,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Temperature (C)") +
  theme_bw() 

# Plot chlorophyll
timeseries_chlorophyll <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = chlorophyll_uL,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Chlorophyll (ug.L)") +
  theme_bw() 

# Plot bacteria
timeseries_bacteria <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = bact,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Est. n of bacterial cells") +
  theme_bw() 

# Plot ammonium
timeseries_ammonium <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = nh4,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Ammonium concentration (umol.L)") +
  theme_bw() 

# Plot nitrite
timeseries_nitrite <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = no2,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Nitrite concentration (umol.L)") +
  theme_bw() 

# Plot nitrate
timeseries_nitrate <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = no3,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Nitrate concentration (umol.L)") +
  theme_bw() 

# Plot phosphate
timeseries_phosphate <- 
  ggplot(time_series, 
       aes(x = visit_id, 
           y = po4,
           shape = vessel,
           colour = subsidy_3,
           group = bromeliad_id)) +
  geom_point(size = 3) +
  geom_line() +
  scale_shape_manual(name = "Vessel",
                     labels = c("Bromeliad", "Wax"),
                     values = c(15, 17)) +
  scale_colour_manual(name = "Subdsidy",
                      labels = c("Litter + Feces", "Litter"),
                      values = c("burlywood", "cadetblue")) + 
  xlab("Hours after subsidy added") +
  ylab("Phosphate concentration (umol.L)") +
  theme_bw() 

# Combine plots
timeseries_grid <- 
  plot_grid(timeseries_temp +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none"),
            timeseries_ph +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none"),
            timeseries_bacteria +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none"),
            timeseries_chlorophyll +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none"),
            timeseries_nitrite +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none"),
            timeseries_nitrate +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none"),
            timeseries_phosphate +
              theme(legend.position = "none"),
            timeseries_ammonium +
              theme(legend.position = "none"),
            ncol = 2)

ggsave(timeseries_grid,
       file = "timeseries.jpg",
       width = 15,
       height = 10)
