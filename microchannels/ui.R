#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(viridis)
library(cowplot)
library(gridGraphics)
library(data.table)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(here)


# Create list of options for buttons in each panel ----------------------------------
# Experiments
exp_choices <- 
  c("Experiment 1", "Experiment 2")

# X axis
x_choices <- 
  c("Weekly measurements", "Time series in exp. 1")

# Facet
facet_choices <- 
  c("Exposition", "Resource", "Bromeliad/Mosquito")

# Y choices
y_choices <- 
  list(HTML(paste0("NO",tags$sub("2"), tags$sup("-"), "(nitrite)")),
       HTML(paste0("NO",tags$sub("3"), tags$sup("-"), "(nitrate)")),
       HTML(paste0("NH",tags$sub("4"), tags$sup("+"), "(ammonium)")),
       "Dissolved inorganic nitrogen (DIN)",
       HTML(paste0("PO",tags$sub("4"), tags$sup("3-"), "(phosphate)")),
       "Bacteria", "Algae", "pH", "Temperature",
       "Mosquito death", "Mosquito pupation", "Mosquito emergence",
       "Time to death", "Time to pupation", "Time to emergence",
       "Larval length at death (mm)", "Dry mass at emergence (mg)",
       "Average wing length of adult (mm)")

# Y values
y_values <- 
  c("NO2", "NO3", "NH4", "DIN", "PO4",
    "Bacteria", "Algae", "pH", "Temperature",
    "Mosquito death", "Mosquito pupation", "Mosquito emergence",
    "Time to death", "Time to pupation", "Time to emergence",
    "Larval length at death (mm)", "Dry mass at emergence (mg)",
    "Average wing length of adult (mm)")


# General dashboard ----------------------------------------------------------------
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
                       ## Experiment
                       radioButtons(
                           inputId = "experiment1",
                           label = "Choose experiment",
                           choices = exp_choices, 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x1",
                           label = "X-axis",
                           choices = x_choices, 
                           inline = TRUE),
                       ## Facet
                       radioButtons(
                           inputId = "facet1",
                           label = "Choose facet",
                           choices = facet_choices, 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                           inputId = "y1",
                           label = "Y-axis",
                           choiceNames = y_choices,
                           choiceValues = y_values,
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
                           choices = exp_choices, 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x2",
                           label = "X-axis",
                           choices = x_choices, 
                           inline = TRUE),
                       ## Facet
                       radioButtons(
                           inputId = "facet2",
                           label = "Choose facet",
                           choices = facet_choices, 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                           inputId = "y2",
                           label = "Y-axis",
                           choiceNames = y_choices,
                           choiceValues = y_values,
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
                           choices = exp_choices, 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x3",
                           label = "X-axis",
                           choices = x_choices, 
                           inline = TRUE),
                       ## Factor 1
                       radioButtons(
                           inputId = "facet3",
                           label = "Choose facet",
                           choices = facet_choices, 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                         inputId = "y3",
                         label = "Y-axis",
                         choiceNames = y_choices,
                         choiceValues = y_values,
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
                           choices = exp_choices, 
                           inline = TRUE),
                       ## X-axis
                       radioButtons(
                           inputId = "x4",
                           label = "X-axis",
                           choices = x_choices, 
                           inline = TRUE),
                       ## Factor 1
                       radioButtons(
                           inputId = "facet4",
                           label = "Choose facet",
                           choices = facet_choices, 
                           inline = TRUE),
                       ## Y-axis
                       radioButtons(
                         inputId = "y4",
                         label = "Y-axis",
                         choiceNames = y_choices,
                         choiceValues = y_values,
                         inline = TRUE)
                   ) # end box 1
            ), # end column 2
        ), # end fluidrow
        

        
        
    ), # end body
    sidebar = dashboardSidebar(disable = TRUE)# here, we only have one tab, so we don't need a sidebar
)
)


