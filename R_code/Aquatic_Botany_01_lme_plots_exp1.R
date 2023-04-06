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

# Extract data from the first three weeks (before)
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





# SessionInfo -------------------------------------------------------------
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 22621)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.1   car_3.1-0       carData_3.0-5   ggeffects_1.1.2 here_1.0.1     
# [6] forcats_0.5.2   stringr_1.4.1   dplyr_1.1.0     purrr_0.3.4     readr_2.1.3    
# [11] tidyr_1.2.0     tibble_3.1.8    ggplot2_3.4.0   tidyverse_1.3.2 lmerTest_3.1-3 
# [16] lme4_1.1-31     Matrix_1.4-1   
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.4.4          jsonlite_1.8.3      splines_4.2.1       modelr_0.1.9       
# [5] assertthat_0.2.1    stats4_4.2.1        tensorA_0.36.2      googlesheets4_1.0.1
# [9] cellranger_1.1.0    pbivnorm_0.6.0      numDeriv_2016.8-1.1 pillar_1.8.1       
# [13] backports_1.4.1     lattice_0.20-45     glue_1.6.2          rvest_1.0.3        
# [17] minqa_1.2.5         colorspace_2.0-3    pkgconfig_2.0.3     broom_1.0.0        
# [21] haven_2.5.0         corpcor_1.6.10      scales_1.2.1        tzdb_0.3.0         
# [25] cubature_2.0.4.5    googledrive_2.0.0   mgcv_1.8-42         generics_0.1.3     
# [29] ellipsis_0.3.2      withr_2.5.0         cli_3.5.0           mnormt_2.1.0       
# [33] crayon_1.5.2        magrittr_2.0.3      readxl_1.4.1        fs_1.5.2           
# [37] fansi_1.0.3         nlme_3.1-157        MASS_7.3-57         xml2_1.3.3         
# [41] tools_4.2.1         hms_1.1.2           gargle_1.2.0        lifecycle_1.0.3    
# [45] reprex_2.0.2        munsell_0.5.0       DHARMa_0.4.5        compiler_4.2.1     
# [49] rlang_1.0.6         grid_4.2.1          nloptr_2.0.3        rstudioapi_0.13    
# [53] lavaan_0.6-12       boot_1.3-28         gtable_0.3.1        abind_1.4-5        
# [57] DBI_1.1.3           R6_2.5.1            lubridate_1.8.0     utf8_1.2.2         
# [61] rprojroot_2.0.3     MCMCglmm_2.33       ape_5.6-2           stringi_1.7.8      
# [65] parallel_4.2.1      Rcpp_1.0.9          vctrs_0.5.2         dbplyr_2.2.1       
# [69] tidyselect_1.2.0    coda_0.19-4   
# 
