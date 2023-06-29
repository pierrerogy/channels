# Load necessary packages
library(tidyverse)
library(ggplot2)
library(ggeffects)
# Nutrient functions written by Kaleigh Davis (https://www.kaleighdavis.com/)

# Nutrients ---------------------------------------------------------------


# Make function to read in plate-shaped csv file, make into dataframe, rename columns
read_plate_function <- function(file_name) {
  data <- read_plate(
    file = file_name,            
    well_ids_column = "Wells")  
  ## PR adds some lines to make sure all sample ids are read as characters
  data <- 
    data %>% 
    mutate_at(names(data[, grepl("_sam", names(data))]), as.character)
}

# Zero each plate according to its own blank values for no2
plate_specific_zeroes_no2 <- function(df){
  standards <- df %>% 
    filter(!is.na(no2_stds))  #isolate no3_standards
  blank_value_no2 <- 
    standards %>% #get blank value
    filter(no2_stds == 0 &
             !is.na(absorbance)) %>% 
    summarise(mean(absorbance)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_no2)
}

# Function to make no3 calibration curves
no2_cal_curve <- function(df) {
  lm(absorbance ~ no2_stds, data = df)
}


# Zero each plate according to its own blank values for no3
plate_specific_zeroes_no3 <- function(df){
  standards <- df %>% 
    filter(!is.na(no3_stds))  #isolate no3_standards
  blank_value_no3 <- standards %>% #get blank value
    filter(no3_stds == 0 &
             !is.na(no3)) %>% 
    summarise(mean(no3)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_no3)
}

# Function to make no3 calibration curves
no3_cal_curve <- function(df) {
  lm(no3 ~ no3_stds, data = df)
}

# Zero each plate according to its own blank values for po4
plate_specific_zeroes_po4 <- function(df){
  standards <- df %>% 
    filter(!is.na(po4_stds))  #isolate no3_standards
  blank_value_po4 <- standards %>% #get blank value
    filter(po4_stds == 0.0 &
             !is.na(absorbance)) %>% 
    summarise(mean(absorbance)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_po4)
}

# Function to make po4 calibration curves
po4_cal_curve <- function(df) {
  lm(absorbance ~ po4_stds, data = df)
}


# Zero each plate according to its own blank values for nh4
plate_specific_zeroes_nh4 <- function(df){
  standards <- df %>% 
    filter(!is.na(nh4_stds))  #isolate standards
  blank_value_nh4 <- standards %>% #get blank value
    filter(nh4_stds == 0.0 &
             !is.na(absorbance)) %>% 
    summarise(mean(absorbance)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_nh4)
}

# Function to make nh44 calibration curves
nh4_cal_curve <- function(df) {
  lm(absorbance ~ nh4_stds, data = df)
}

# Zero each plate according to its own blank values for bacteria
plate_specific_zeroes_bact <- function(df){
  standards <- df %>% 
    filter(!is.na(bact_stds))  #isolate standards
  blank_value_bact <- standards %>% #get blank value
    filter(bact_stds == 0.0 &
             !is.na(absorbance)) %>% 
    summarise(mean(absorbance)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_bact)
}

# Function to make bact4 calibration curves
bact_cal_curve <- function(df) {
  lm(absorbance ~ bact_stds, data = df)
}



# Models ------------------------------------------------------------------
# Get axis label
axis_label <- function(parameter, axis){
  # Many forks to get what we want
  if(stringr::str_detect(string = parameter, pattern = "din"))
    ## x
    if(axis == "x")
    {ret <- 
      expression("DIN concentration (log"~mu*"mol"*".L"^"-1"*",scaled and centered)")}
    ## y
  else if(axis == "y"){
    ret <- 
      expression("DIN concentration (log"~mu*"mol"*".L"^"-1"*")")}

  if(stringr::str_detect(string = parameter, pattern = "po4"))
    ## x
    if(axis == "x")
    {ret <- 
      expression("PO"["4"]^"3-"*" concentration (log"~mu*"mol"*".L"^"-1"*", scaled and centered)")}
    ## y
  else if(axis == "y"){
    ret <- 
      expression("PO"["4"]^"3-"*" concentration (log"~mu*"mol"*".L"^"-1"*")")}

  if(parameter == "np_log")
    ## x
    if(axis == "x")
    {ret <- 
      "N:P ratio (log,s scaled and centered)"}
  ## y
  else if(axis == "y"){
    ret <- 
      "N:P ratio (log)"}
  
  
  if(parameter == "pH")
    ret <- 
      "pH"
  
  if(parameter == "temperature_C")
    ret <- 
      "Temperature (°C),"
  
  if(stringr::str_detect(string = parameter, pattern = "chloro"))
    ## x
    if(axis == "x")
    {ret <- 
      expression("Chlorophyll-a concentration (log"~mu*"g"*".L"^"-1"*", scaled and centered)")}
  ## y
  else if(axis == "y"){
    ret <- 
      expression("Chlorophyll-a concentration (log"~mu*"g"*".L"^"-1"*")")}

  
  if(stringr::str_detect(string = parameter, pattern = "bact"))
    ## x
    if(axis == "x")
    {ret <- 
      expression("Bacteria concentration (log x"*" 10"^"12"*""*".L"^"-1"*", scaled and centered)")}
  ## y
  else if(axis == "y"){
    ret <- 
      expression("Bacteria concentration (log x"*" 10"^"12"*""*".L"^"-1"*")")}
  
  if(stringr::str_detect(string = parameter, pattern = "exposure"))
    ret <- 
      "Light exposure"
  
  if(stringr::str_detect(string = parameter, pattern = "moz"))
    ret <- 
      "Number of larvae \nin mesocosm (log)"
  
  if(stringr::str_detect(string = parameter, pattern = "mass"))
    ret <- 
      "Dry adult biomass (log mg)"
  
  
  if(stringr::str_detect(string = parameter, pattern = "death"))
    ret <- 
      "Age at death (log days)"
  
  if(parameter == "time_death")
    ret <- 
      "Age at death (log days)"
  
  
  if(parameter == "time_pupation")
    ret <- 
      "Age at pupation (log days)"
  
  
  if(parameter == "time_emergence")
    ret <- 
      "Age at emergence (log days)"
  
  if(parameter == "death")
    ret <- 
      "Individual survival success"
  
  
  if(parameter == "pup")
    ret <- 
      "Individual pupation success"
  
  
  if(parameter == "emergence")
    ret <- 
      "Individual emergence success"
  
  
  if(stringr::str_detect(string = parameter, pattern = "wing"))
    ret <- 
      "Average wing length (log mm)"
  
  
  if(stringr::str_detect(string = parameter, pattern = "size"))
    ret <- 
      "Body length at death (log mm)"
  
  # Return
  return(ret)
}

# Make plots for models 
plot_model_nice <- function(model, xax, yax, scale = "none", type, data){
  # Prepare data
    if(yax %in% c("bact_log", "chlorophyll_ugL_log", "n_moz")){
      model_effect <- 
        purrr::flatten_df(brms::conditional_effects(model,
                                                    effects = c(xax),
                                                    method = "fitted")[1])
    } else 
      {
        model_effect <- 
          brms::conditional_effects(model,
                                    effects = c("subsidy:exposure"),
                                    method = "fitted")$`subsidy:exposure`
        }
  
  # Change n_moz name
  if(yax == "n_moz")
  {yax = "n_moz_log"}
  
  # Get axes label
  ylab <- 
    axis_label(parameter = yax,
               axis = "y")
  xlab <- 
    axis_label(parameter = xax,
               axis = "x")
  
  # Plot according to model type
  if(type == "points"){
    ret <- 
      ggplot(data = model_effect,
             aes(x = exposure, 
                 y = estimate__), 
             colour = subsidy) + 
      geom_point(size = 3,
                 aes(colour = subsidy),
                 position = position_dodge(0.5)) +
      geom_errorbar(aes(ymin = lower__, 
                        ymax = upper__,
                        colour = subsidy), 
                    width = 0.2,
                    position = position_dodge(0.5)) +
      ggtitle("") +
      xlab(xlab) +
      scale_x_discrete(labels = c("Exposed", "Shaded")) +
      ylab(ylab) +
      scale_color_manual(name = "Subsidy",
                         labels = c("Litter", "Litter + feces"), 
                         values = c("tan1", "tan4")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    ## Modify scale if logged
    if(scale == "log"){
      ret <- 
        ret +
        geom_jitter(data = data,
                    mapping = aes(x = exposure,
                                  y = log(get(yax)),
                                  colour = subsidy),
                    alpha = 0.3)
    } else if(scale == "prob")
    {ret <- 
      ret +
      geom_point(aes(x = exposure, 
                     y = get(yax),
                     colour = subsidy), 
                 data = data, 
                 position = position_jitter(w = 0.8, h = 0),
                 alpha = 0.3) 
    } else
      {ret <- 
        ret +
        geom_jitter(data = data,
                    mapping = aes(x = exposure,
                                  y = get(yax),
                                  colour = subsidy),
                    alpha = 0.3) 
      }
  }
  
  ## Another fork depending on which group of microorganism
  if(type == "lines"){
    ## Make mini loop to frame depending on exposed or shaded microcosms
   ret <- 
     ggplot(data = model_effect,
                   aes(x = effect1__, 
                       y = estimate__)) + 
            geom_line(colour = "burlywood3") +
            geom_ribbon(aes(ymin = lower__, 
                            ymax = upper__,
                            fill = "burlywood3"), 
                        colour = NA,
                        alpha = 0.2) +
            geom_jitter(size = 3,
                       data = data,
                       aes(colour = subsidy,
                           x = get(xax),
                           y = get(yax)),
                       height = 0.5) +
            guides(fill="none") +
            ggtitle("") +
            xlab(xlab) +
            ylab(ylab) +
            scale_color_manual(name = "Subsidy",
                               labels = c("Litter", "Litter + feces"), 
                               values = c("tan1", "tan4")) +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black")) 
  }
 
    # Return
    return(ret)
}

# Get equivalent models from loo output
get_eq_models <- function(loo_output){
  ret <- 
    ## Extract diff table values
    tibble::tibble(model = row.names(loo_output$diffs),
                   elpd_diff = loo_output$diffs[,1],
                   se_diff = loo_output$diffs[,2]) %>% 
    ## Get ratio of elpd/se to get equivalent models
    dplyr::mutate(ratio = abs(elpd_diff/se_diff)) %>% 
    ## Keep best and equivalent models
    dplyr::filter(ratio < 2 | is.nan(ratio))
  
  # Return
  return(ret)
}