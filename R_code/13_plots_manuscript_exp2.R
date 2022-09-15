# Path analysis
# Load libraries
library(tidyverse)
library(here)
library(ggplot2)
library(cowplot)
source("R_code/functions.R")


# Figure 2 - Models with effects of nutrients on all microorganisms ----------------------------------------------------------------
## Get coefficients for algae
chloro_exp2_bactnut_coefs <- 
  get_coefs(data = exp2_center,
            model = chloro_exp2_bactnut,
            bacteria = T)

## DIN - algae
{
  ## Plot
  figure2a <- 
    ggplot(data = chloro_exp2_bactnut_coefs,
           aes(x = din_scale,
               y = fit,
               colour = shading,
               fill = shading)) +
    geom_smooth(method = 'glm') +
    geom_jitter(data = exp2_center,
                aes(x = din_scale,
                    y = chlorophyll_ugL,
                    colour = shading)) +
    scale_y_continuous(trans = "log",
                       breaks = c(1, 20, 50, 400)) +
    xlab(expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))) +
    ylab(expression(paste("Chlorophyll-a concentration (",mu, "g."*"L"^"-1"*")"))) +
    scale_colour_manual(name = "Light exposure",
                        labels = c("Exposed", "Shaded"), 
                        values = c("goldenrod", "grey50")) +
    scale_fill_manual(name = "Light exposure",
                      labels = c("Exposed", "Shaded"), 
                      values = c("goldenrod", "grey50")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  }

## Bacteria - algae
{
  ## Plot
  figure2b <- 
    ggplot(data = chloro_exp2_bactnut_coefs,
           aes(x = bact_scale,
               y = fit,
               colour = shading,
               fill = shading)) +
    geom_smooth(method = 'glm') +
    geom_jitter(data = exp2_center,
                aes(x = bact_scale,
                    y = chlorophyll_ugL,
                    colour = shading)) +
    scale_y_continuous(trans = "log",
                       breaks = c(1, 20, 50, 400)) +
    xlab(expression("Bacteria concentration (x"*"10"^"12"*""*".L"^"-1"*")")) +
    ylab(expression(paste("Chlorophyll-a concentration (", mu ,"g"*".L"^"-1"*")"))) +
    scale_colour_manual(name = "Light exposure",
                        labels = c("Exposed", "Shaded"), 
                        values = c("goldenrod", "grey50")) +
    scale_fill_manual(name = "Light exposure",
                      labels = c("Exposed", "Shaded"), 
                      values = c("goldenrod", "grey50")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  
  
}

## Get legend
legend_2 <- 
  cowplot::get_legend(figure2a)

# DIN - bacteria
{## Get coefficients
  bact_exp2_all_coefs <- 
    get_coefs(data = exp2_center,
              model = bact_exp2_all,
              chloro = T)
  
  ## Plot
  figure2c <- 
    ggplot(data = bact_exp2_all_coefs,
           aes(x = din_scale,
               y = fit,
               colour = shading,
               fill = shading)) +
    geom_smooth(method = 'glm') +
    geom_jitter(data = exp2_center,
                aes(x = din_scale,
                    y = bact,
                    colour = shading)) +
    scale_y_continuous(trans = "log",
                       breaks = c(0, 1, 6)) +
    xlab(expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))) +
    ylab(expression("Bacteria concentration  (x"*"10"^"12"*""*".L"^"-1"*")")) +
    scale_colour_manual(name = "Light exposure",
                        labels = c("Exposed", "Shaded"), 
                        values = c("goldenrod", "grey50")) +
    scale_fill_manual(name = "Light exposure",
                      labels = c("Exposed", "Shaded"), 
                      values = c("goldenrod", "grey50")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  ## Plot
  figure2d <- 
    ggplot(data = bact_exp2_all_coefs,
           aes(x = chloro_scale,
               y = fit,
               colour = shading,
               fill = shading)) +
    geom_smooth(method = 'glm') +
    scale_y_continuous(trans = "log",
                       breaks = c(1, 6, 30))+
    geom_jitter(data = exp2_center,
                aes(x = chloro_scale,
                    y = bact,
                    colour = shading)) +
    xlab(expression(paste("Chlorophyll-a concentration (", mu, "g"*".L"^"-1"*")"))) +
    ylab(expression("Bacteria concentration  (x"*"10"^"12"*""*".L"^"-1"*")")) +
    scale_colour_manual(name = "Light exposure",
                        labels = c("Exposed", "Shaded"), 
                        values = c("goldenrod", "grey50")) +
    scale_fill_manual(name = "Light exposure",
                      labels = c("Exposed", "Shaded"), 
                      values = c("goldenrod", "grey50")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  }

# Treatment - bacteria
{## Get coefficients
  bact_exp2_all_coefs <- 
    get_coefs(data = exp2_center,
              model = bact_exp2_all,
              feces = T)
  
  ## Plot
  figure2e <- 
    ggplot(data = bact_exp2_all_coefs,
           aes(x = subsidy,
               y = fit,
               colour = shading)) +
    geom_point(position = position_dodge(width = 0.5),
               size = 3) +
    geom_errorbar(aes(ymin = lwr,
                      ymax = upr,
                      colour = shading),
                  position = position_dodge(width = 0.5)) +
    geom_jitter(data = exp2_center,
                aes(x = subsidy,
                    y = bact,
                    colour = shading),
                alpha = 0.3) +
    scale_x_discrete(name = "Terrestrial subsidy",
                     labels = c("Litter only", "Litter and feces")) +
    ylab(expression("Bacteria concentration  (x"*"10"^"12"*""*".L"^"-1"*")")) +
    scale_colour_manual(name = "Light exposure",
                        labels = c("Exposed", "Shaded"), 
                        values = c("goldenrod", "grey50")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
}


## Combine plots
{
  figure2 <- 
    cowplot::plot_grid(figure2a +
                         ggtitle("a") +
                         theme(legend.position = "none"),
                       figure2b +
                         ggtitle("b") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       legend_2,
                       figure2c +
                         ggtitle("c") +
                         theme(legend.position = "none"),
                       figure2d +
                         ggtitle("d") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       figure2e +
                         ggtitle("e") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       ncol = 3,
                       nrow = 2,
                       rel_widths = c(1, 0.8, 0.8))
}

## Save plot
ggplot2::ggsave(figure2,
                height = 7,
                width = 10,
                filename = paste0(here::here("figures",
                                             "exp2_figure2.jpeg")))


# Figure 3 - Models with effects of nutrients on mosquitoes ----------------------------------------------------------------
# DIN - mosquitoes
{### Effect
  moznut_din <- 
    as.data.frame(ggeffects::ggpredict(moznut_model,
                                       terms = c("din_scale"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    #dplyr::filter(group == "exposed") %>% 
    dplyr::select(-group)
  ### Plot
  figure3a <- 
    ggplot2::ggplot(data = moznut_din,
                    aes(x = x,
                        y = predicted)) +
    geom_line(colour = "lightcyan3",
              size = 2) +
    geom_ribbon(aes(ymin = conf.low, 
                    ymax = conf.high, 
                    alpha = 0.2),
                colour = NA,
                fill = "lightcyan3") +
    scale_y_continuous(trans = "log10",
                       breaks = c(0.0001, 1, 6),
                       labels = c("0", "1", "6")) +
    geom_jitter(data = exp2_center %>% 
                  dplyr::mutate(n_moz = n_moz + 0.0001),
                aes(x = din_scale,
                    y = n_moz),
                colour = "lightcyan3") +
    xlab(expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))) +
    ylab("Number of live larvae in cup") +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

# PO4 - mosquitoes
{### Effect
  moznuteffect_po4 <- 
    as.data.frame(ggeffects::ggpredict(moznut_model,
                                       terms = c("po4_scale"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::select(-group)
  ### Plot
  figure3b <- 
    ggplot2::ggplot(data = moznuteffect_po4,
                    aes(x = x,
                        y = predicted)) +
    geom_line(colour = "lightcyan3",
              size = 2) +
    geom_ribbon(aes(ymin = conf.low, 
                    ymax = conf.high, 
                    alpha = 0.2),
                colour = NA,
                fill = "lightcyan3") +
    scale_y_continuous(trans = "log10",
                       breaks = c(0.0001, 1, 6),
                       labels = c("0", "1", "6")) +
    geom_jitter(data = exp2_center %>% 
                  dplyr::mutate(n_moz = n_moz + 0.0001),
                aes(x = po4_scale,
                    y = n_moz),
                colour = "lightcyan3") +
    xlab(expression(paste("PO"["4"]^"3-"*" concentration (", mu, "mol"*".L"^"-1"*")"))) +
    ylab("Number of live larvae in cup") +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
}

## Combine plots
{
  figure3 <- 
    cowplot::plot_grid(figure3a +
                         ggtitle("a") +
                         theme(legend.position = "none"),
                       figure3b +
                         ggtitle("b") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       ncol = 2,
                       nrow = 1)
}

## Save plot
ggplot2::ggsave(figure3,
                height = 4,
                width = 8,
                filename = paste0(here::here("figures",
                                             "exp2_figure3.jpeg")))

# Figure 4 - Mosquito survival growth and health -----------------------------------
## Combine plots
figure4 <- 
    cowplot::plot_grid(avg_timedeath_plot +
                         ggtitle("a") +
                         theme(legend.position = "none"),
                       avg_sizedeath_plot +
                         ggtitle("b") +
                         theme(legend.position = "none"),
                       avg_winglength_plot  +
                         ggtitle("c") +
                         theme(legend.position = "none"),
                       tot_biomass_plot +
                         ggtitle("d") +
                         theme(legend.position = "none"),
                       ncol = 2,
                       nrow = 2)

## Get legend
legend_4 <- 
  cowplot::get_legend(sizedeath_plot)

## Add legend to grid
figure4 <- 
  cowplot::plot_grid(figure4, 
                     legend_4, 
                     ncol = 2, 
                     rel_widths = c(1, 0.2))

## Save plot
ggplot2::ggsave(figure4,
                height = 6,
                width = 8,
                filename = paste0(here::here("figures",
                                             "exp2_figure4.jpeg")))



# Figure S1 - Treatments on nutrients -------------------------------------
## Combine plots
figures1 <- 
  cowplot::plot_grid(exp2_din_reslight_plot +
                       ggtitle("a") +
                       theme(legend.position = "none"),
                     exp2_po4_reslight_plot+
                       ggtitle("b") +
                       theme(legend.position = "none"),
                     exp2_np_reslight_plot +
                       ggtitle("c") +
                       theme(legend.position = "none"),
                     exp2_pH_reslight_plot+
                       ggtitle("d") +
                       theme(legend.position = "none"),
                     exp2_temperature_reslight_plot+
                       ggtitle("e") +
                       theme(legend.position = "none"),
                     nrow = 3,
                     ncol = 2,
                     rel_heights = c(1, 1, 1))

## Add legend to grid
figures1 <- 
  cowplot::plot_grid(figures1, 
                     legend_4, 
                     ncol = 2, 
                     rel_widths = c(1, 0.2))

## Save plot
ggplot2::ggsave(figures1,
                height = 8,
                width = 8,
                filename = paste0(here::here("figures",
                                             "exp2_figures1.jpeg")))

# Figure S2 - Other mosquito survival, growth and health (indiv) --------
## Combine plots
figures2 <- 
  cowplot::plot_grid(avg_timepupation_plot +
                       ggtitle("a") +
                       theme(legend.position = "none"),
                     avg_timeemergence_plot +
                       ggtitle("b") +
                       theme(legend.position = "none"),
                     prop_probdeath_plot +
                       ggtitle("c") +
                       theme(legend.position = "none"),
                     prop_probpup_plot +
                       ggtitle("d") +
                       theme(legend.position = "none"),
                     prop_probemergence_plot  +
                       ggtitle("e") +
                       theme(legend.position = "none"),
                     ncol = 2,
                     nrow = 3)


## Add legend to grid
figures2 <- 
  cowplot::plot_grid(figures2, 
                     legend_4, 
                     ncol = 2, 
                     rel_widths = c(1, 0.2))

## Save plot
ggplot2::ggsave(figures2,
                height = 10,
                width = 8,
                filename = paste0(here::here("figures",
                                             "exp2_figures2.jpeg")))
