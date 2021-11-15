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
  

# Define server logic required to draw a plots
shinyServer(function(input, output) {
    
    # Read in the prepared data 
    exp1 <- 
        read.csv(here::here("appdata", 
                            "bromeliad_tax_exp.csv")) %>% 
        dplyr::rename(temperature_C = temp,
                      bromsquito = vessel,
                      pH = ph)
    
    exp2 <- 
        read.csv(here::here("appdata", 
                            "weekly_measurements_exp2.csv")) %>% 
      dplyr::rename(bromsquito = larvae)
    
    mosquitoes <- 
        read.csv(here::here("appdata", 
                            "mosquitoes.csv"))
    
    
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
                facet_par = input$facet3,
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
          if(input$y1 %notin% c("Mosquito pupation", "Mosquito death")){
            lineplot1 <-
              line_blank_plot(lineplot1,  plot1_dats())} else
                if(input$y1 %in% c("Mosquito pupation", "Mosquito death")){
                  lineplot1 <-
                    point_blank_plot(lineplot1,  plot1_dats())}

            
        }
        
        
        
        # Print the plot with correct labels
        lineplot1 +
          ylab(get_y_label(input$y1)) +
          xlab(get_x_label(input$x1))
        
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
        if(input$y2 %notin% c("Mosquito pupation", "Mosquito death")){
          lineplot2 <-
            line_blank_plot(lineplot2,  plot2_dats())} else
              if(input$y2 %in% c("Mosquito pupation", "Mosquito death")){
                lineplot2 <-
                  point_blank_plot(lineplot2,  plot2_dats())}
        
        
        
        
      }
      
      
      
      # Print the plot with correct labels
      lineplot2 +
        ylab(get_y_label(input$y2)) +
        xlab(get_x_label(input$x2))
        
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
        if(input$y3 %notin% c("Mosquito pupation", "Mosquito death")){
          lineplot3 <-
            line_blank_plot(lineplot3,  plot3_dats())} else
              if(input$y3 %in% c("Mosquito pupation", "Mosquito death")){
                lineplot3 <-
                  point_blank_plot(lineplot3,  plot3_dats())}

        
      }
      
      
      
      # Print the plot with correct labels
      lineplot3 +
        ylab(get_y_label(input$y3)) +
        xlab(get_x_label(input$x3))
        
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
        if(input$y4 %notin% c("Mosquito pupation", "Mosquito death")){
          lineplot4 <-
            line_blank_plot(lineplot4,  plot4_dats())} else
              if(input$y4 %in% c("Mosquito pupation", "Mosquito death")){
                lineplot4 <-
                  point_blank_plot(lineplot4,  plot4_dats())}

        
      }
      
      
      
      # Print the plot with correct labels
      lineplot4  +
        ylab(get_y_label(input$y4)) +
        xlab(get_x_label(input$x4))
        
    }) # end renderplot command
    
    
    
})
