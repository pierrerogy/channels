# Path analysis
# Load libraries
library(tidyverse)
library(here)
library(ggplot2)
library(cowplot)
source("R_code/functions.R")


# Figure 2 - Models with effects of nutrients on all microorganisms ----------------------------------------------------------------
# Get legend
legend_2 <- 
  cowplot::get_legend(plot1)

# Organise plots
figure2 <- 
    cowplot::plot_grid(plot45[[1]] +
                         ggtitle("a") +
                         theme(legend.position = "none"),
                       plot45[[2]] +
                         ggtitle("b") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       plot67[[1]] +
                         ggtitle("c") +
                         theme(legend.position = "none"),
                       plot67[[2]]  +
                         ggtitle("d") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       plot89[[1]] +
                         ggtitle("e") +
                         theme(legend.position = "none"),
                       plot89[[2]]  +
                         ggtitle("f") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       plot1011[[1]] +
                         ggtitle("g") +
                         theme(legend.position = "none"),
                       plot1011[[2]]  +
                         ggtitle("h") +
                         theme(legend.position = "none",
                               axis.title.y = element_blank()),
                       ncol = 2)

# Add legend
figure2 <- 
  cowplot::plot_grid(figure2,
                     legend_2,
                     ncol = 2,
                     rel_widths = c(0.9, 0.1))

# Save plot
ggplot2::ggsave(figure2,
                height = 11,
                width = 10,
                bg = "white",
                filename = paste0(here::here("figures",
                                             "exp2_figure2.jpeg")))


# Figure 3 - Mosquito survival growth and health -----------------------------------
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
