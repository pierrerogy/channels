# Path analys for experiment 1

# Load libraries
library(tidyverse)
library(here)
library(lavaan)
library(semPlot)
library(lme4)
library(ggeffects)
library(car)
library(cowplot)

# Load data
# Overall data
exp1 <- 
  read.csv(here::here("microchannels",
                      "appdata", 
                      "bromeliad_tax_exp.csv")) %>% 
  dplyr::rename(temperature_C = temp,
                pH = ph) %>% 
  ## Remove time series
  dplyr::mutate(visit_id = as.numeric(visit_id),
                bromeliad_id = as.numeric(bromeliad_id)) %>% 
  dplyr::filter(!is.na(visit_id)) %>% 
  ## Make visit_character
  dplyr::mutate(visit_id = as.character(visit_id)) %>% 
  ## Sum nitrogen values into dissolved inorganic nitrogen
  dplyr::mutate(DIN = no2 + no3 + nh4) %>% 
  ## Put capital letter for vessel columns
  dplyr::mutate(vessel = stringr::str_to_title(vessel))

# Extract data from before
## Before
exp1_before <- 
  exp1 %>% 
  ## Keep relevant visit ids
  dplyr::filter(visit_id %in% c("1", "2", "3")) %>% 
  ## Remove columns we will not use
  dplyr::select(-pH, -temperature_C, -larvae, -subsidy_id,
                -capacity_mL, -longest_leaf_length_mm,
                -n_leaves, -subsidy, -no2, -no3, -nh4, -date) %>% 
  dplyr::mutate(subsidy_1 = ifelse(subsidy_1 == "litter_feces",
                                   "Litter and feces",
                                   "Litter only"))
## Make data wide format for lavaan
exp1_before_wide <- 
  exp1_before %>% 
  ## Remove all NAs
  dplyr::filter(!is.na(po4) & !is.na(DIN) & 
                  !is.na(bact) & !is.na(chlorophyll_ugL)) %>% 
  ## Join scaled data
  dplyr::left_join(exp1_before %>%
                     ## Scale and rename data
                     dplyr::mutate_at(vars(po4, DIN, bact, chlorophyll_ugL),
                                      ~as.numeric(scale(.))) %>% ## make numeric to remove attributes
                     dplyr::rename_at(vars(po4, DIN, bact, chlorophyll_ugL),
                                      ~ paste( .x, "scale",sep = "_")),
                   by = c("visit_id", "bromeliad_id",
                          "vessel", "shading", "subsidy_1")) %>% 
  ## Pivot wider
  tidyr::pivot_wider(names_from = visit_id,
                     values_from = c(DIN, po4, chlorophyll_ugL, bact,
                                     DIN_scale, po4_scale, chlorophyll_ugL_scale, bact_scale)) %>%
  ## Add code for subsidy
  dplyr::mutate(subsidy_1 = ifelse(subsidy_1 == "Litter only",
                                 0, 1),
                vessel = ifelse(vessel == "Wax",
                                1, 0),
                shading = ifelse(shading == "light",
                                 0, 1))

# Path analysis ----------------------------------------------------------
# Model description 
SEM_before <- '
    # Structural relation
      chlorophyll_ugL_scale_2 ~ b1*shading + b2*subsidy_1 + b3*DIN_scale_2 + b4*po4_scale_2+b31*vessel
      bact_scale_2 ~ b5*shading + b6*vessel + b7*subsidy_1 + b8*DIN_scale_2 + b9*po4_scale_2
      DIN_scale_2 ~ b10*vessel + b11*subsidy_1
      po4_scale_2 ~ b12*vessel + b13*subsidy_1 
      chlorophyll_ugL_scale_3 ~ b14*shading + b15*subsidy_1 + b16*DIN_scale_3 + b17*po4_scale_3 + b18*chlorophyll_ugL_scale_2+ b32*vessel
      bact_scale_3 ~ b19*shading + b20*vessel + b21*subsidy_1 + b22*DIN_scale_3 + b23*po4_scale_3 + b24*bact_scale_2
      DIN_scale_3 ~ b25*vessel + b26*subsidy_1 + b27*DIN_scale_2
      po4_scale_3 ~ b28*vessel + b29*subsidy_1 + b30*po4_scale_2
    # Variance structure of exogenous variables
      vessel ~~ vessel
      subsidy_1 ~~ subsidy_1
      shading ~~ shading
    # Residual covariance
      po4_scale_2 ~~ DIN_scale_2
      po4_scale_3 ~~ DIN_scale_3
      chlorophyll_ugL_scale_2 ~~ bact_scale_2
      chlorophyll_ugL_scale_3 ~~ bact_scale_3
    '

# Fit model    
fit_SEM_before <-
  lavaan::sem(SEM_before,
              data = exp1_before_wide,
              estimator = "ML",
              meanstructure = F)

# # Explore possible extra paths to improve fit. All fit measures not great. Added a few paths in model above. Not more paths left plausible.
# View(lavaan::modindices(fit_SEM_before,
#                         standardized = TRUE,
#                         sort. = TRUE))

# Summarize output, include standardized estimates and fit measures  
summary(fit_SEM_before, 
        rsquare = TRUE, 
        standardized = TRUE, 
        fit.measures = TRUE)
# Here we want:
# user p value above 0.05
# CFI above 0.95 
# TLI above 0.9
# RMSEA under 0.08
# SRMR under 0.08
#

# Visualize 
semPlot::semPaths(fit_SEM_before, 
                 shapeMan="rectangle",sizeMan = 14,sizeMan2 = 5,
                 rotation = 2,layout="tree2",what="std",
                 posCol="black",edge.width=0.5,
                 style="Lisrel",edge.label.position=0.3,
                 fade = F,
                 edge.label.cex=1)



# Linear models -----------------------------------------------------------
# DIN
hist(log(exp1_before$DIN))
# Model
din_model_before <- 
  lme4::lmer(log(DIN) ~ shading*subsidy_1*vessel + 
               (1|visit_id) + (1|bromeliad_id),
             data = exp1_before)
# Assumptions
plot(din_model_before)
# Tests
car::Anova(din_model_before)
# Plot
## Data
din_effect_before <- 
  ggeffects::ggpredict(din_model_before,
                       terms = c("shading", "vessel", "subsidy_1"),
                       type = "re",
                       ci.level = 0.95)
## Plot
din_plot_before <- 
  plot(din_effect_before,
       ci = T,
       dot.size = 3,
       line.size = 1,
       rawdata = T) + 
  ggtitle("") + 
  ylab(expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))) +
  xlab("Light exposure") +
  scale_colour_manual(name = "Vessel",
                      labels = c("Bromeliad", "Wax"), 
                      values = c("darkgreen", "gray")) +
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

# PO4
hist(log(exp1_before$po4))
# Model
po4_model_before <- 
  lme4::lmer(log(po4) ~ shading*subsidy_1*vessel + 
               (1|visit_id) + (1|bromeliad_id),
             data = exp1_before)
# Assumptions
plot(po4_model_before)
# Tests
car::Anova(po4_model_before)
# Plot
## Data
po4_effect_before <- 
  ggeffects::ggpredict(po4_model_before,
                       terms = c("shading", "vessel", "subsidy_1"),
                       type = "re",
                       ci.level = 0.95)
## Plot
po4_plot_before <- 
  plot(po4_effect_before,
       ci = T,
       dot.size = 3,
       line.size = 1,
       rawdata = T) + 
  ggtitle("") + 
  ylab(expression(paste("PO"["4"]^"3-"*" concentration (", mu, "mol"*".L"^"-1"*")"))) +
  xlab("Light exposure") +
  scale_colour_manual(name = "Vessel",
                      labels = c("Bromeliad", "Wax"), 
                      values = c("darkgreen", "gray")) +
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

# Bacteria
hist(log(exp1_before$bact))
# Model
bact_model_before <- 
  lme4::lmer(log(bact) ~ shading*subsidy_1*vessel + 
               (1|visit_id) + (1|bromeliad_id),
             data = exp1_before)
# Assumptions
plot(bact_model_before)
# Tests
car::Anova(bact_model_before)
# Plot
## Data
bact_effect_before <- 
  ggeffects::ggpredict(bact_model_before,
                       terms = c("shading", "vessel", "subsidy_1"),
                       type = "re",
                       ci.level = 0.95)
## Plot
bact_plot_before <- 
  plot(bact_effect_before,
       ci = T,
       dot.size = 3,
       line.size = 1,
       rawdata = T) + 
  ggtitle("") + 
  ylab(expression("Number of bacterial cells (x"*"10"^"12"*""*".L"^"-1"*")")) +
  xlab("Light exposure") +
  scale_colour_manual(name = "Vessel",
                      labels = c("Bromeliad", "Wax"), 
                      values = c("darkgreen", "gray")) +
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

# Chlorophyll
hist(log(exp1_before$chlorophyll_ugL))
# Model
chloro_model_before <- 
  lme4::lmer(log(chlorophyll_ugL) ~ shading*subsidy_1*vessel + 
               (1|visit_id) + (1|bromeliad_id),
             data = exp1_before)
# Assumptions
plot(chloro_model_before)
# Tests
car::Anova(chloro_model_before)
# Plot
## Data
chloro_effect_before <- 
  ggeffects::ggpredict(chloro_model_before,
                       terms = c("shading", "vessel", "subsidy_1"),
                       ci.level = 0.95)
## Plot
chloro_plot_before <- 
  plot(chloro_effect_before,
       ci = T,
       dot.size = 3,
       line.size = 1,
       rawdata = T) + 
  ggtitle("") + 
  scale_y_continuous(trans = "log",
                     breaks = c(2, 20, 80, 140)) +
  ylab(expression(paste("Chlorophyll-a concentration (", mu ,"g"*".L"^"-1"*")"))) +
  xlab("Light exposure") +
  scale_colour_manual(name = "Vessel",
                      labels = c("Bromeliad", "Wax"), 
                      values = c("darkgreen", "gray")) +
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


# Get legends
legend <- 
  cowplot::get_legend(bact_plot_before)
# Combine figures
figure2 <- 
  cowplot::plot_grid(din_plot_before +
                       ggtitle("A") +
                       theme(legend.position = "none") +
                       xlab(""),
                     po4_plot_before +
                       ggtitle("B") +
                       theme(legend.position = "none") +
                       xlab(""),
                     chloro_plot_before +
                       ggtitle("C") +
                       theme(legend.position = "none"),
                      bact_plot_before +
                       ggtitle("D") +
                       theme(legend.position = "none"),
                     ncol = 2,
                     nrow = 2)

# Combine everything
figure_2 <- 
  cowplot::plot_grid(figure2, 
                     legend, 
                     ncol = 2,
                     rel_widths = c(1, 0.2))
# Save figure
ggsave(here::here("figures",
                  "exp1_figure2.jpg"),
       figure_2,
       width = 10,
       height=10)



