# Functions for app

# Libraries
library(dplyr)
library(stringr)

# Select data -------------------------------------------------------------
get_those_dats <- function(y, x, facet_par, experiment, 
                           exp1, exp2, mosquitoes){
  # Make facet selection flexible
  facet_var <- 
    ifelse(facet_par == "Exposition",
           "shading", ifelse(facet_par == "Resource",
                      "subsidy", "bromsquito"))
  col_var <- 
    ifelse(facet_par == "Exposition",
           "bromsquito", ifelse(facet_par == "Resource",
                             "shading", "subsidy"))
  
  lty_var <- 
    ifelse(facet_par == "Exposition",
           "subsidy", ifelse(facet_par == "Resource",
                             "bromsquito", "shading"))
  
  # Select correct experiment
  if(experiment == "Experiment 1"){
     ## Here is where the x can either be time lapse or across weeks
    if(x == "Time series in exp. 1"){
      dats <- 
        exp1 %>% 
        dplyr::filter(stringr::str_detect(visit_id, "H")) %>% 
        dplyr::rename(x = visit_id,
                      sample_id = bromeliad_id,
                      facet = eval(facet_var),
                      col = eval(col_var),
                      lty = eval(lty_var)) %>% 
        dplyr::mutate(x = stringr::str_replace_all(x, "H", ""),
                      x = as.numeric(x))
    
    }
    if(x == "Weekly measurements"){
      dats <- 
        exp1 %>% 
        dplyr::filter(!stringr::str_detect(visit_id, "H")) %>% 
        dplyr::rename(x = date,
                      sample_id = bromeliad_id,
                      facet = eval(facet_var),
                      col = eval(col_var),
                      lty = eval(lty_var))
      
    }
    
  }
  
  if(experiment == "Experiment 2"){
    dats <- 
      exp2 %>% 
      dplyr::rename(x = date,
                    sample_id = cup_number,
                    facet = eval(facet_var),
                    col = eval(col_var),
                    lty = eval(lty_var))
    
    
  }
  # If ph is asked
  if(y == "pH"){
    dats <- 
      dats %>% 
      dplyr::rename(y = pH)
    
  } else
    
  # If no2 is asked
  if(y == "Temperature"){
    dats <- 
      dats %>% 
      dplyr::rename(y = temperature_C)
      
    } else
  
  
  # If no2 is asked
  if(y == "NO2 (nitrite)"){
    dats <- 
      dats %>% 
      dplyr::rename(y = no2)
      
  } else
  
  # If no3 are asked
  if(y == "NO3 (nitrate)"){
    dats <- 
      dats %>% 
      dplyr::rename(y = no3)
  } else
  
  # If nh4 are asked
  if(y == "NH4 (ammonium)"){
    dats <- 
      dats %>% 
      dplyr::rename(y = nh4)
    
  } else
  
  # If po4 is asked
  if(y == "PO4 (phosphate)"){
    dats <- 
      dats %>% 
      dplyr::rename(y = po4)
  } else
  
  # If bacterial biomass is asked
  if(y == "Bacteria"){
    dats <- 
      dats %>% 
      dplyr::rename(y = bact)
  } else
  
  # If chlorophyll is asked
  if(y == "Algae"){
      dats <- 
        dats %>% 
        dplyr::rename(y = chlorophyll_ugL)
    } else
  
  # If mosquito death is asked
  if(y == "Mosquito death"){
        dats <- 
          mosquitoes %>% 
          dplyr::filter(!is.na(death)) %>% 
          dplyr::select(cup_number, death, shading, subsidy, larvae) %>% 
          dplyr::group_by(death, cup_number, shading, subsidy, larvae) %>% 
          dplyr::tally() %>% 
          dplyr::rename(x = death,
                        y = n,
                        sample_id = cup_number)
      } else  
        
  # If mosquito pupation is asked
  if(y == "Mosquito pupation"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(pupation)) %>% 
      dplyr::select(cup_number, pupation, shading, subsidy, larvae) %>% 
      dplyr::group_by(pupation, cup_number, shading, subsidy, larvae) %>% 
      dplyr::tally() %>% 
      dplyr::rename(x = pupation,
                    y = n,
                    sample_id = cup_number)
    
  }  else
  
  # If mosquito time to death
  if(y == "Time to death"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(time_death)) %>% 
      dplyr::select(cup_number, time_death, shading, subsidy, larvae) %>% 
      dplyr::group_by(time_death, cup_number, shading, subsidy, larvae) %>% 
      dplyr::tally() %>% 
      dplyr::rename(x = time_death,
                    y = n) %>% 
      ## Add blank column
      dplyr::mutate(sample_id = NA)
    
  } else
  
  # If mosquito time to emergence
  if(y == "Time to emergence"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(time_emergence)) %>% 
      dplyr::select(cup_number, time_emergence, shading, subsidy, larvae) %>% 
      dplyr::group_by(time_emergence, cup_number, shading, subsidy, larvae) %>% 
      dplyr::tally() %>% 
      dplyr::rename(x = time_emergence,
                    y = n) %>% 
      ## Add blank column
      dplyr::mutate(sample_id = NA)
    
  } 
  
  # Return data
  return(dats) 
  
}



# Select y axis label -----------------------------------------------------
get_y_label <- function(y){
  # If ph is asked
  if(y == "pH"){
    lab <- 
      "pH"
    
  } else
    
    # If t is asked
    if(y == "Temperature"){
      lab <- 
        "Temperature (C)"
      
    } else
      
      
      # If no2 is asked
      if(y == "NO2 (nitrite)"){
        lab <- 
          "NO2 (umol/L)"
        
      } else
        
        # If no3 is asked
        if(y == "NO3 (nitrate)"){
          lab <- 
            "NO3 (umol/L)"
        } else
          
          # If nh4 is asked
          if(y == "NH4 (ammonium)"){
            lab <- 
              "NH4(umol/L)"
          } else
            
            # If po4 is asked
            if(y == "PO4 (phosphate)"){
              lab <- 
                "PO4 (umol/L)"
            } else
              
              # If bacterial biomass is asked
              if(y == "Bacteria"){
                lab <- 
                  "Approx. n. of bacterial cells (x 10e12/L)"
              } else
                
                # If chlorophyll is asked
                if(y == "Algae"){
                  lab <- 
                    "Chlorophyll (ug.L)"
                } else
                  
                  # If mosquito death is asked
                  if(y == "Mosquito death"){
                    lab <- 
                      "n dead larvae found per cup"
                  } else  
                    
                    # If mosquito pupation is asked
                    if(y == "Mosquito pupation"){
                      lab <- 
                        "n pupae found per cup"
                      
                    }  else
                      
                      # If time to death is asked
                      if(y == "Time to death"){
                        lab <- 
                          "n"
                        
                      }  else
                        
                        # If time to emergence is asked
                        if(y == "Time to emergence"){
                          lab <- 
                            "n"
                          
                        }

  
  
  # Return label
  return(lab) 
  
  
  
  
  
}




# Select x axis label -----------------------------------------------------
get_x_label <- function(x, y){
  
  # If weekly measurements are asked
  if(x == "Weekly measurements"){
    lab <-
      "Day of the year"
  }
  # Small exception with time to death/emergence
  if(stringr::str_detect(y, "Time")){
    lab <-
      "Number of days"
  }
  
  # If time series is asked
  if(x == "Time series in exp. 1"){
    lab <-
      "Hours after subsidy added"
  }
      
  
  # Return label
  return(lab) 
  
  
  
  
  
}




# Make blank plot ---------------------------------------------------------
blank_plot <- function(dats){
  ## Plot with nothing but 
  return(ggplot(data = dats,
         aes(x = x,
             y = y,
             group = sample_id)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())) 
  
  
}


# Add line to  blank plot -------------------------------------------------
line_blank_plot <- function(blank_plot, dats, y){
  
  ## Colour, linetype and facet all depend on the data
  plot <- 
    blank_plot +
    geom_line(data = dats,
              aes(x = x,
                  y = y,
                  color = col,
                  linetype = lty)) +
    facet_wrap(~ facet) +
    ## Remove automatic labels on legend
    scale_colour_discrete(name = NULL) +
    scale_linetype_discrete(name = NULL)
  
  ## Now some plots need to be logged, depending on y
  if(y %in% c("NH4 (ammonium)", "Algae")){
    plot <- 
      plot +
      scale_y_continuous(trans = "log10")
    
  }
    
  
  
  
  return(plot)
  
  
    
  
  
}

# Add points to  blank plot -------------------------------------------------
point_blank_plot <- function(blank_plot, dats){
  ## Colour and shape fixed
  return(blank_plot +
           geom_point(data = dats,
                      position = position_dodge(width = 0.9),
                      aes(x = x,
                          y = jitter(y),
                          color = shading,
                          shape = subsidy),
                      size = 2) +
           ylim(0, 6) +
           ## Remove automatic labels on legend
           scale_colour_discrete(name = NULL) +
           scale_shape_discrete(name = NULL))
  
  
}
