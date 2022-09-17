# Path analys for experiment 1

# Load libraries
library(tidyverse)
library(here)
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
car::Anova(din_model_before,
           type = "3")
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
  scale_colour_manual(name = "Waxing",
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
car::Anova(po4_model_before,
           type = "3")
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
  scale_colour_manual(name = "Waxing",
                      labels = c("Unwaxed", "Waxed"), 
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
car::Anova(bact_model_before,
           type = "3")
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
  scale_colour_manual(name = "Waxing",
                      labels = c("Unwaxed", "Waxed"), 
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
car::Anova(chloro_model_before,
           type ="3")
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
    scale_colour_manual(name = "Waxing",
                      labels = c("Unwaxed", "Waxed"), 
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



# Combine plots -----------------------------------------------------------
# Get legends
legend <- 
  cowplot::get_legend(bact_plot_before)
# Combine figures
figure1 <- 
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
figure_1 <- 
  cowplot::plot_grid(figure1, 
                     legend, 
                     ncol = 2,
                     rel_widths = c(1, 0.2))
# Save figure
ggsave(here::here("figures",
                  "exp1_figure1.jpg"),
       figure_1,
       width = 10,
       height=10,
       bg = "white") ## somehow came out with black background



