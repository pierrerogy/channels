# Functions for app

# Libraries
library(dplyr)
library(stringr)


# Select between facet, colour and line parameters ------------------------
select_plot_pars <- function(facet_par, experiment){
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
  
  # Give facet names 
  if(experiment == "Experiment 1")
  {if(facet_var == "shading")
    facet_labs <- c(light = "Light", 
                    shade = "Shade") else
      if(facet_var == "subsidy")
        facet_labs <- c(litter_feces = "Litter and feces",
                        litter_only = "Litter only") else
          if(facet_var == "bromsquito" & experiment == "Experiment 1") 
            facet_labs <- c(bromeliad = "Bromeliad",
                            wax = "Wax")}
  
  if(experiment == "Experiment 2")
  {if(facet_var == "shading")
    facet_labs <- c(exposed= "Light", 
                    shaded = "Shade") else
                      if(facet_var == "subsidy")
                        facet_labs <- c(litter_feces = "Litter and feces",
                                        litter = "Litter only") else
                                          if(facet_var == "bromsquito") 
                                            facet_labs <- c(absent = "Larvae absent",
                                                            present = "Larvae present")}
  
                              
                              
                              
  
  # Give col parameters
  if(col_var == "shading")
    col_params <- c("Light", "Shade") else
      if(col_var == "subsidy")
        col_params <-  c("Litter and feces", "Litter only") else
          if(col_var == "bromsquito" & experiment == "Experiment 1") 
            col_params <- c("Bromeliad", "Wax") else
              col_params <- c("Larvae absent", "Larvae present")
  
  # Give col values
  if(col_var == "shading")
    col_vals <- c("goldenrod", "grey50") else
      if(col_var == "subsidy")
        col_vals <-  c("tan4", "tan1") else
          if(col_var == "bromsquito" & experiment == "Experiment 1") 
            col_vals <- c("darkgreen", "gray") else
              col_vals <- c("gray", "black")
  
  # Give lty parameters
  if(lty_var == "shading")
    lty_params <- c("Light", "Shade") else
      if(lty_var == "subsidy")
        lty_params <-  c("Litter and feces", "Litter only") else
          if(lty_var == "bromsquito" & experiment == "Experiment 1") 
            lty_params <- c("Bromeliad", "Wax") else
              lty_params <- c("Larvae absent", "Larvae present")
  
  
  # Return values
  return(list(facet_var, col_var, lty_var, facet_labs,
           col_params, col_vals, lty_params))
  
}


# Select data -------------------------------------------------------------
get_those_dats <- function(y, x, facet_par, experiment, 
                           exp1, exp2, mosquitoes){
  
  # Get plot parameters
  plot_pars <- 
    select_plot_pars(facet_par, experiment)
  
  # Select correct experiment
  if(experiment == "Experiment 1"){
     ## Here is where the x can either be time lapse or across weeks
    if(x == "Time series in exp. 1"){
      dats <- 
        exp1 %>% 
        dplyr::filter(stringr::str_detect(visit_id, "H")) %>% 
        dplyr::rename(x = visit_id,
                      sample_id = bromeliad_id,
                      facet = eval(plot_pars[[1]]),
                      col = eval(plot_pars[[2]]),
                      lty = eval(plot_pars[[3]])) %>% 
        dplyr::mutate(x = stringr::str_replace_all(x, "H", ""),
                      x = as.numeric(x))
    
    }
    if(x == "Weekly measurements"){
      dats <- 
        exp1 %>% 
        dplyr::filter(!stringr::str_detect(visit_id, "H")) %>% 
        dplyr::rename(x = date,
                      sample_id = bromeliad_id,
                      facet = eval(plot_pars[[1]]),
                      col = eval(plot_pars[[2]]),
                      lty = eval(plot_pars[[3]]))
      
    }
    
  }
  
  if(experiment == "Experiment 2"){
    dats <- 
      exp2 %>% 
      dplyr::rename(x = date,
                    sample_id = cup_number,
                    facet = eval(plot_pars[[1]]),
                    col = eval(plot_pars[[2]]),
                    lty = eval(plot_pars[[3]]))
    
    
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
  if(y == "NO2"){
    dats <- 
      dats %>% 
      dplyr::rename(y = no2)
      
  } else
  
  # If no3 is asked
  if(y == "NO3"){
    dats <- 
      dats %>% 
      dplyr::rename(y = no3)
  } else
  
  # If nh4 is asked
  if(y == "NH4"){
    dats <- 
      dats %>% 
      dplyr::rename(y = nh4)
    
  } else
    
  # If DIN is asked
  if(y == "DIN"){
    dats <- 
      dats %>% 
      dplyr::rename(y = DIN)
    
  } else
  
  # If po4 is asked
  if(y == "PO4"){
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
  if(y == "Mosquito death" & experiment == "Experiment 2"){
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
  if(y == "Mosquito pupation" & experiment == "Experiment 2"){
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
  
  # If mosquito emergence is asked
  if(y == "Mosquito emergence" & experiment == "Experiment 2"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(emergence)) %>% 
      dplyr::select(cup_number, emergence, shading, subsidy, larvae) %>% 
      dplyr::group_by(emergence, cup_number, shading, subsidy, larvae) %>% 
      dplyr::tally() %>% 
      dplyr::rename(x = emergence,
                    y = n,
                    sample_id = cup_number)
    
  }  else  
    
  # If mosquito time to death
  if(y == "Time to death" & experiment == "Experiment 2"){
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
    
  # If mosquito time to pupation
  if(y == "Time to pupation" & experiment == "Experiment 2"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(time_pupation)) %>% 
      dplyr::select(cup_number, time_pupation, shading, subsidy, larvae) %>% 
      dplyr::group_by(time_pupation, cup_number, shading, subsidy, larvae) %>% 
      dplyr::tally() %>% 
      dplyr::rename(x = time_pupation,
                    y = n) %>% 
      ## Add blank column
      dplyr::mutate(sample_id = NA)
    
  } else  
  
  # If mosquito time to emergence
  if(y == "Time to emergence" & experiment == "Experiment 2"){
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
  
  # If size of larvae at death
  if(y == "Larval length at death (mm)" & experiment == "Experiment 2"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(size_mm)) %>% 
      dplyr::select(cup_number, time_death, size_mm, shading, subsidy, larvae) %>% 
      dplyr::rename(x = time_death,
                    y = size_mm) %>% 
      ## Add blank column
      dplyr::mutate(sample_id = NA)
    
  } 
  
  # If dry mass of emerging adults
  if(y == "Dry mass at emergence (mg)" & experiment == "Experiment 2"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(dry_mass_mg)) %>% 
      dplyr::select(cup_number, time_emergence, dry_mass_mg, shading, subsidy, larvae) %>% 
      dplyr::rename(x = time_emergence,
                    y = dry_mass_mg) %>% 
      ## Add blank column
      dplyr::mutate(sample_id = NA)
    
  } 
  
  # If wing length of emerged adults
  if(y == "Average wing length of adult (mm)" & experiment == "Experiment 2"){
    dats <- 
      mosquitoes %>%
      dplyr::filter(!is.na(wing_length)) %>% 
      dplyr::select(cup_number, time_emergence, wing_length, shading, subsidy, larvae) %>% 
      dplyr::rename(x = time_emergence,
                    y = wing_length) %>% 
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
      if(y == "NO2"){
        lab <- 
          expression(paste("NO"["2"]^"-"*" concentration (", mu, "mol"*".L"^"-1"*")"))
        
      } else
        
        # If no3 is asked
        if(y == "NO3"){
          lab <- 
            expression(paste("NO"["3"]^"-"*" concentration (", mu, "mol"*".L"^"-1"*")"))
        } else
          
          # If nh4 is asked
          if(y == "NH4"){
            lab <- 
              expression(paste("NH"["4"]^"+"*" concentration (", mu, "mol"*".L"^"-1"*")"))
          } else
            
            # If DIN is asked
            if(y == "DIN"){
              lab <- 
                expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))
            } else
            
            # If po4 is asked
            if(y == "PO4"){
              lab <- 
                expression(paste("PO"["4"]^"3-"*" concentration (", mu, "mol"*".L"^"-1"*")"))
            } else
              
              # If bacterial biomass is asked
              if(y == "Bacteria"){
                lab <- 
                  expression("Number of bacterial cells (x"*"10"^"12"*""*".L"^"-1"*")")
              } else
                
                # If chlorophyll is asked
                if(y == "Algae"){
                  lab <- 
                    expression(paste("Chlorophyll-a concentration (", mu ,"g"*".L"^"-1"*")"))
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
                      
                      # If mosquito emegence is asked
                      if(y == "Mosquito emergence"){
                        lab <- 
                          "n emerged adults per cup"
                        
                      }  else
                      
                      # If time to death/pupation/emergence is asked
                      if(stringr::str_detect(y, "Time")){
                        lab <- 
                          "n"
                        
                      }  else
                        
                        # If size of larvae at death
                        if(y == "Larval length at death (mm)"){
                            lab <- 
                              "Individual larval length at death (mm)"
                          
                        } else
                          
                          # If dry mass of emerging adults
                          if(y == "Dry mass at emergence (mg)"){
                            lab <- 
                              "Individual dry mass at emergence (mg)"
                            
                          } else
                            
                            # If wing length of emerged adults
                            if(y == "Average wing length of adult (mm)"){
                              lab <- 
                                "Individual average wing length of adult (mm)"
                              
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
  # Exception with time to death/emergence
  if(stringr::str_detect(y, "Time")){
    lab <-
      paste("Number of days until", 
             stringr::str_remove(y, "Time to "))
  }
  
  # Exception for size at death
  if(y == "Larval length at death (mm)"){
    lab <-
      "Number of days until death"
  }
  
  # Exception for dry mass and wing length
  if(y %in% c("Dry mass at emergence (mg)", "Average wing length of adult (mm)")){
    lab <-
      "Number of days until emergence"
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

# Add line to blank plot -------------------------------------------------
line_blank_plot <- function(blank_plot, dats, y, facet_par, experiment){
  
  # Get plotting parameters
  plot_pars <- 
    select_plot_pars(facet_par, experiment)
  
  ## Colour, linetype and facet all depend on the data
  plot <- 
    blank_plot +
    geom_line(data = dats,
              aes(x = x,
                  y = y,
                  color = col,
                  linetype = lty)) +
    facet_wrap(~ facet,
               labeller = as_labeller(plot_pars[[4]])) +
    ## Remove automatic labels on legend
    scale_colour_manual(name = NULL,
                        labels = plot_pars[[5]],
                        values = plot_pars[[6]]) +
    scale_linetype_manual(name = NULL,
                          labels = plot_pars[[7]],
                          values = c(1,2))
  
  ## Now some plots need to be logged, depending on y
  if(y %in% c("NH4", "Algae")){
    plot <- 
      plot +
      scale_y_continuous(trans = "log10")
    
  }
    
  
  
  
  return(plot)
  
  
    
  
  
}

# Add points to blank plot -------------------------------------------------
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
           ylim(0, NA) +
           ## Remove automatic labels on legend
           scale_colour_manual(name = NULL,
                                 labels = c("Exposed", "Shaded"),
                                 values = c("goldenrod", "grey50")) +
           scale_shape_manual(name = NULL,
                              labels = c("Litter only", "Litter and feces"),
                              values = c(17, 16)))
  
  
}
