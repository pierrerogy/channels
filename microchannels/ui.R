#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
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


# Define UI for application that draws plots of the data
shinyUI(dashboardPage(
    skin = "purple",
    header=
        dashboardHeader(
            title = "Microchannels experiments",
            titleWidth = 500 # since we have a long title, we need to extend width element in pixels
        ),
    
    # create dashboard body - this is the major UI element
    body= dashboardBody(
        
        # First row has inputs for plots 1 and 2 ----------------------------------------------------
        fluidRow(
            # 1/2 of page (width = 6 of 12 columns)
            column(width = 6,
                   box(width = NULL, status = "primary", collapsible = T,
                       title  = "Inputs for plot 1", solidHeader = T,
                       ## Change colour of box headers
                       tags$style(HTML(
                       ".box.box-solid.box-primary>.box-header {
                       color:#fff;
                       background:#555299;
                       background-color: #555299;
                       }")),
                       ## Landscape
                       radioButtons(
                           inputId = "experiment1",
                           label = "Choose experiment",
                           choices = c("Experiment 1", "Experiment 2"), 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x1",
                           label = "X-axis",
                           choices = c("Weekly measurements", "Time series in exp. 1"), 
                           inline = TRUE),
                       ## Facet
                       radioButtons(
                           inputId = "facet1",
                           label = "Choose facet",
                           choices = c("Exposition", "Resource", "Bromeliad/Mosquito"), 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                           inputId = "y1",
                           label = "Y-axis",
                           choices = c("NO2 (nitrite)", "NO3 (nitrate)",
                                       "NH4 (ammonium)", "PO4 (phosphate)",
                                       "Bacteria", "Algae", "pH", "Temperature",
                                       "Mosquito death", "Mosquito pupation",
                                       "Time to death", "Time to emergence"
                                       ),
                           inline = TRUE)
                   ) # end box 1
            ), # end column 1
            # 1/2 of page (width = 6 of 12 columns)
            column(width = 6,
                   box(width = NULL, status = "primary", collapsible = T,
                       title  = "Inputs for plot 2", solidHeader = T,
                       ## Landscape
                       radioButtons(
                           inputId = "experiment2",
                           label = "Choose experiment",
                           choices = c("Experiment 1", "Experiment 2"), 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x2",
                           label = "X-axis",
                           choices = c("Weekly measurements", "Time series in exp. 1"), 
                           inline = TRUE),
                       ## Facet
                       radioButtons(
                           inputId = "facet2",
                           label = "Choose facet",
                           choices = c("Exposition", "Resource", "Bromeliad/Mosquito"), 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                           inputId = "y2",
                           label = "Y-axis",
                           choices = c("NO2 (nitrite)", "NO3 (nitrate)",
                                       "NH4 (ammonium)", "PO4 (phosphate)",
                                       "Bacteria", "Algae", "pH", "Temperature",
                                       "Mosquito death", "Mosquito pupation",
                                       "Time to death", "Time to emergence"),
                           inline = TRUE)
                   ) # end box 1
            ) # end column 2
            # end column 4
        ), # end fluidrow
        # Second row has plots 1 and 2 ----------------------------------------------------
        fluidRow(column(width = 6,
                        box(width = NULL, status = "primary",
                            solidHeader = TRUE, 
                            title = "Plot1",
                            plotOutput("plot1", 
                                       height = 350))),
                 column(width = 6,
                        box(width = NULL, status = "primary",
                            solidHeader = TRUE, 
                            title = "Plot2",
                            plotOutput("plot2", 
                                       height = 350)))
                 
        ),
        
        # Third row has plots 3 and 4 ----------------------------------------------------
        fluidRow(column(width = 6,
                        box(width = NULL, status = "primary",
                            solidHeader = TRUE, 
                            title = "Plot3",
                            plotOutput("plot3", 
                                       height = 350))),
                 column(width = 6,
                        box(width = NULL, status = "primary",
                            solidHeader = TRUE, 
                            title = "Plot4",
                            plotOutput("plot4", 
                                       height = 350)))
                 
        ),
        
        # Fourth row has inputs for plots 3 and 4 ---------------------------------
        fluidRow(
            # 1/2 of page (width = 6 of 12 columns)
            column(width = 6,
                   box(width = NULL, status = "primary", collapsible = T,
                       title  = "Inputs for plot 3", solidHeader = T,
                       ## Landscape
                       radioButtons(
                           inputId = "experiment3",
                           label = "Choose experiment",
                           choices = c("Experiment 1", "Experiment 2"), 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x3",
                           label = "X-axis",
                           choices = c("Weekly measurements", "Time series in exp. 1"), 
                           inline = TRUE),
                       ## Factor 1
                       radioButtons(
                           inputId = "facet3",
                           label = "Choose facet",
                           choices = c("Exposition", "Resource", "Bromeliad/Mosquito"), 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                           inputId = "y3",
                           label = "Y-axis",
                           choices = c("NO2 (nitrite)", "NO3 (nitrate)",
                                       "NH4 (ammonium)", "PO4 (phosphate)",
                                       "Bacteria", "Algae", "pH", "Temperature",
                                       "Mosquito death", "Mosquito pupation",
                                       "Time to death", "Time to emergence"),
                           inline = TRUE)
                   ) # end box 1
            ), # end column 1
            # 1/2 of page (width = 6 of 12 columns)
            column(width = 6,
                   box(width = NULL, status = "primary", collapsible = T,
                       title  = "Inputs for plot 4", solidHeader = T,
                       ## Landscape
                       radioButtons(
                           inputId = "experiment4",
                           label = "Choose experiment",
                           choices = c("Experiment 1", "Experiment 2"), 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x4",
                           label = "X-axis",
                           choices = c("Weekly measurements", "Time series in exp. 1"), 
                           inline = TRUE),
                       ## Factor 1
                       radioButtons(
                           inputId = "facet4",
                           label = "Choose facet",
                           choices = c("Exposition", "Resource", "Bromeliad/Mosquito"), 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                           inputId = "y4",
                           label = "Y-axis",
                           choices = c("NO2 (nitrite)", "NO3 (nitrate)",
                                       "NH4 (ammonium)", "PO4 (phosphate)",
                                       "Bacteria", "Algae", "pH", "Temperature",
                                       "Mosquito death", "Mosquito pupation",
                                       "Time to death", "Time to emergence"),
                           inline = TRUE)
                   ) # end box 1
            ), # end column 2
        ), # end fluidrow
        

        
        
    ), # end body
    sidebar = dashboardSidebar(disable = TRUE)# here, we only have one tab, so we don't need a sidebar
)
)

