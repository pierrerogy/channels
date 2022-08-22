#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(tidyr)
library(gridGraphics)
library(stringr)
library(data.table)
library(shiny)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(here)
source(here::here("functions.R"))
`%notin%` <-Negate(`%in%`)
  

# Lists of mosquito variables unique to experiment 2 -------------------------------
mosq_vars <- 
  c("Mosquito death", "Mosquito pupation", "Mosquito emergence",
    "Time to death", "Time to pupation", "Time to emergence",
    "Larval length at death (mm)", "Dry mass at emergence (mg)",
    "Average wing length of adult (mm)")


# Build server ------------------------------------------------------------
# Define server logic required to draw a plots
shinyServer(function(input, output) {
    
    # Read in the prepared data 
  exp1 <- 
        read.csv(here::here("appdata", 
                            "bromeliad_tax_exp.csv")) %>% 
        dplyr::rename(temperature_C = temp,
                      bromsquito = vessel,
                      pH = ph) %>% 
        ## Make dissolved inorganic nitrogen column
        dplyr::mutate(DIN = no2 + no3 + nh4) %>% 
        ## Remove bromeliads in second part of the experiment
        ## that received the wrong subsidy
        dplyr::filter(visit_id %in% c("1", "2", "3", "H0", "H1",
                                      "H4", "H8", "H24", "H48") |
                        (visit_id %in% c("4", "5", "6") & 
                           subsidy %in% c("litter_only_litter_only",
                                          "litter_feces_litter_feces"))) %>% 
        ## Put value of subsidy in subsidy_1 for time series data
        dplyr::mutate(subsidy_1 = ifelse(is.na(subsidy_1),
                                         subsidy, subsidy_1)) %>% 
        ## Change name of subsidy column
        dplyr::select(-subsidy) %>% 
        dplyr::rename(subsidy = subsidy_1)
    
    exp2 <- 
        read.csv(here::here("appdata", 
                            "weekly_measurements_exp2.csv")) %>% 
      dplyr::rename(bromsquito = larvae) %>% 
      ## Make dissolved inorganic nitrogen column, make subisdy category same than exp1
      dplyr::mutate(DIN = no2 + no3 + nh4,
                    subsidy = ifelse(subsidy == "litter",
                                     "litter_only", subsidy))
    
    mosquitoes <- 
        read.csv(here::here("appdata", 
                            "mosquitoes.csv")) %>% 
      ## Add average of two wings
      dplyr::mutate(wing_length = (left_wing_mm + right_wing_mm)/2)
    
    
    # Make reactive datasets ---------------------------------------------------
    # Subset data for each plot depending on selection
    plot1_dats <- reactive({
        
        ## Get data
        dats <- 
            get_those_dats(
                y = input$y1, 
                x = input$x1, 
                facet_par = input$facet1, 
                experiment = input$experiment1,
                exp1 = exp1,
                exp2 = exp2, 
                mosquitoes = mosquitoes)
        
        ## Return data
        return(dats)
        
    })
    plot2_dats <- reactive({
        
        ## Get data
        dats <- 
            get_those_dats(
                y = input$y2, 
                x = input$x2, 
                facet_par = input$facet2,
                experiment = input$experiment2,
                exp1 = exp1,
                exp2 = exp2, 
                mosquitoes = mosquitoes)
        
        ## Return data
        return(dats)
        
    })
    plot3_dats <- reactive({
        
        ## Get data
        dats <- 
            get_those_dats(
                y = input$y3, 
                x = input$x3, 
                facet_par = input$facet3,
                experiment = input$experiment3,
                exp1 = exp1,
                exp2 = exp2, 
                mosquitoes = mosquitoes)
        
        ## Return data
        return(dats)
    })
    plot4_dats <- reactive({
        
        ## Get data
        dats <- 
            get_those_dats(
                y = input$y4, 
                x = input$x4, 
                facet_par = input$facet4,
                experiment = input$experiment4,
                exp1 = exp1,
                exp2 = exp2, 
                mosquitoes = mosquitoes)
        
        ## Return data
        return(dats)
    })
    
    
    # Make plots --------------------------------------------------------------
    
    output$plot1 <- renderPlot({
      # Make blank plot  
      lineplot1 <-
            blank_plot(plot1_dats()) +
            ggtitle(input$experiment1)
        
        # IF a data object exists, update the blank ggplot.
        # basically this makes it not mess up when nothing is selected
        
        if(nrow(plot1_dats()) > 1){
          ## Plot type depends on y
          if(input$y1 %notin% mosq_vars){
            lineplot1 <-
              line_blank_plot(lineplot1,  plot1_dats(), input$y1, input$facet1, input$experiment1)} else
                if(input$y1 %in% mosq_vars){
                  lineplot1 <-
                    point_blank_plot(lineplot1,  plot1_dats())}

            
        }
        
        
        
        # Print the plot with correct labels
        lineplot1 +
          ylab(get_y_label(input$y1)) +
          xlab(get_x_label(input$x1, input$y1))
        
    }) # end renderplot command
    
    output$plot2 <- renderPlot({
      # Make blank plot   
      lineplot2 <-
        blank_plot(plot2_dats()) +
        ggtitle(input$experiment2)
      
      # IF a data object exists, update the blank ggplot.
      # basically this makes it not mess up when nothing is selected
      
      if(nrow(plot2_dats()) > 1){
        ## Plot type depends on y
        if(input$y2 %notin% mosq_vars){
          lineplot2 <-
            line_blank_plot(lineplot2,  plot2_dats(), input$y2, input$facet2, input$experiment2)} else
              if(input$y2 %in% mosq_vars){
                lineplot2 <-
                  point_blank_plot(lineplot2,  plot2_dats())}
        
        
        
        
      }
      
      
      
      # Print the plot with correct labels
      lineplot2 +
        ylab(get_y_label(input$y2)) +
        xlab(get_x_label(input$x2, input$y2))
        
    }) # end renderplot command
    
    output$plot3 <- renderPlot({
      # Make blank plot  
      lineplot3 <-
        blank_plot(plot3_dats()) +
        ggtitle(input$experiment3)
      
      # IF a data object exists, update the blank ggplot.
      # basically this makes it not mess up when nothing is selected
      
      if(nrow(plot3_dats()) > 1){
        ## Plot type depends on y
        if(input$y3 %notin% mosq_vars){
          lineplot3 <-
            line_blank_plot(lineplot3,  plot3_dats(), input$y3, input$facet3, input$experiment3)} else
              if(input$y3 %in% mosq_vars){
                lineplot3 <-
                  point_blank_plot(lineplot3,  plot3_dats())}

        
      }
      
      
      
      # Print the plot with correct labels
      lineplot3 +
        ylab(get_y_label(input$y3)) +
        xlab(get_x_label(input$x3, input$y3))
        
    }) # end renderplot command
    
    output$plot4 <- renderPlot({
      # Make blank plot  
      lineplot4 <-
        blank_plot(plot4_dats()) +
        ggtitle(input$experiment4)
      
      # IF a data object exists, update the blank ggplot.
      # basically this makes it not mess up when nothing is selected
      
      if(nrow(plot4_dats()) > 1){
        ## Plot type depends on y
        if(input$y4 %notin% mosq_vars){
          lineplot4 <-
            line_blank_plot(lineplot4,  plot4_dats(), input$y4, input$facet4, input$experiment4)} else
              if(input$y4 %in% mosq_vars){
                lineplot4 <-
                  point_blank_plot(lineplot4,  plot4_dats())}

        
      }
      
      
      
      # Print the plot with correct labels
      lineplot4  +
        ylab(get_y_label(input$y4)) +
        xlab(get_x_label(input$x4, input$y4))
        
    }) # end renderplot command
    
    
    
})
