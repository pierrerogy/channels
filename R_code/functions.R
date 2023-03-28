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
# Get y axis label
y_axis_label <- function(parameter){
  # Many forks to get what we want
  if(parameter == "din")
    ret <- 
      expression(paste("DIN concentration (", mu, "mol"*".L"^"-1"*")"))
  if(parameter == "po4")
    ret <- 
      expression(paste("PO"["4"]^"3-"*" concentration (", mu, "mol"*".L"^"-1"*")"))
  if(parameter == "np")
    ret <- 
      "N:P ratio"
  if(parameter == "pH")
    ret <- 
      "pH"
  if(parameter == "temperature_C")
    ret <- 
      "Temperature (°C)"
  if(parameter %in% c("chlorophyll_ugL", "bact_log_scale"))
    ret <- 
      expression(paste("Chlorophyll-a concentration (", mu ,"g"*".L"^"-1"*")"))
  if(parameter %in% c("po4_log_scale", "din_log_scale", "chlorophyll_ugL_log_scale"))
    ret <- 
      expression("Bacteria concentration  (x"*"10"^"12"*""*".L"^"-1"*")")
  
  # Return
  return(ret)
}

# Get x axis label
x_axis_label <- function(parameter){
  # Three forks
  if(!stringr::str_detect(string = parameter, pattern = "scale"))
    ret <- 
      "Light exposure"
  if(parameter == "din_log_scale")
    ret <- expression(paste("Scaled DIN concentration (", mu, "mol"*".L"^"-1"*")"))
  if(parameter == "po4_log_scale")
    ret <- 
      expression(paste("Scaled PO"["4"]^"3-"*" concentration (", mu, "mol"*".L"^"-1"*")"))
  if(parameter == "bact_log_scale")
    ret <- 
      expression("Scaled bacteria concentration  (x"*"10"^"12"*""*".L"^"-1"*")")
  if(parameter == "chlorophyll_ugL_log_scale")
    ret <- 
      expression(paste("Scaled chlorophyll-a concentration (", mu ,"g"*".L"^"-1"*")"))
 
  # Return as expression
  return(ret)
}

# Generate mock data to predict values for complex models
generate_mock_data <- function(parameter){
  # First, generate data over the range of values of the parameter
  ret <- 
   tibble::tibble(parameter = rep(seq(min(eval(parse(text = paste0("exp2_center$", parameter)))), 
                                      max(eval(parse(text = paste0("exp2_center$", parameter)))),
                                      by = 0.1),
                                  times = 4))
  
  # Add custom number of values of treatment values, and rando values for others
  ret <- 
    ret %>% 
    dplyr::mutate(exposure = as.factor(rep(c("exposed", "shaded"),
                                           each = nrow(ret) / 2)),
                  subsidy = as.factor(rep(rep(c("litter", "litter_feces"),
                                              each = nrow(ret) / 4),
                                          times = 2))) %>% 
    dplyr::bind_cols(bact_log = 1:nrow(ret),
                     week = 1:nrow(ret),
                     cup_number = 1:nrow(ret))
  
  # Three way fork depending on the independent variable
  if(parameter == "din_log_scale"){
    ret <- 
      ret %>% 
      dplyr::rename(din_log_scale = parameter) %>% 
      dplyr::mutate(chlorophyll_ugL_log_scale = mean(exp2_center$chlorophyll_ugL_log_scale),
                    po4_log_scale = mean(exp2_center$po4_log_scale))} else
                    if(parameter == "chlorophyll_ugL_log_scale"){
                      ret <- 
                        ret %>% 
                        dplyr::rename(chlorophyll_ugL_log_scale = parameter) %>% 
                        dplyr::mutate(din_log_scale = mean(exp2_center$din_log_scale),
                                      po4_log_scale = mean(exp2_center$po4_log_scale))} else
                                      {ret <- 
                                        ret %>% 
                                        dplyr::rename(po4_log_scale = parameter) %>% 
                                        dplyr::mutate(chlorophyll_ugL_log_scale = mean(exp2_center$chlorophyll_ugL_log_scale),
                                                      din_log_scale = mean(exp2_center$din_log_scale))}

  # Get predicted values
  ret <- 
    ret %>% 
    dplyr::bind_cols(predict.MCMCglmm(model, 
                                      newdata = ret,
                                      marginal = model$Random$formula,
                                      type = "response",
                                      level = 0.95,
                                      interval = "prediction"))  %>% 
    dplyr::rename(predicted = fit,
                  conf.low = lwr,
                  conf.high = upr)  
  
  # Return predicted values
  return(ret)
  
}

# Make plots for models 
plot_model_nice <- function(model, parameter, scale = "none", type){
  # Prepare data
  if(parameter == "bact_log_scale"){
    model_effect <- 
      tibble::as_tibble(ggeffects::ggpredict(model,
                                         terms = c("bact_log_scale", "exposure", "subsidy"),
                                         type = "re",
                                         ci.level = 0.95)) %>% 
      dplyr::rename(bact_log_scale = x,
                    exposure = group,
                    subsidy = facet)
    
  } else
    if(parameter %in% c("chlorophyll_ugL_log_scale", "po4_log_scale", "din_log_scale")){
      model_effect <- 
        generate_mock_data(parameter = parameter)
    } else 
      {
  model_effect <- 
    as.data.frame(ggeffects::ggpredict(model = model,
                                       terms = c("exposure", "subsidy"),
                                       type = "re",
                                       ci.level = 0.95)) %>% 
    dplyr::rename(exposure = x,
                  subsidy = group)}
  
  # Transform data according to scale
  if(scale == "log")
    model_effect <- 
      model_effect %>% 
      dplyr::mutate(across(c(predicted, conf.low, conf.high),
                           ~ exp(.)))
  if(scale == "sqrt")
    model_effect <- 
      model_effect %>% 
      dplyr::mutate(across(c(predicted, conf.low, conf.high),
                           ~ .x^2))
  
  # Get axes label
  ylab <- 
    y_axis_label(parameter)
  xlab <- 
    x_axis_label(parameter)
  
  # Plot according to model type
  if(type == "points"){
    ret <- 
      ggplot(data = model_effect,
             aes(x = exposure, 
                 y = predicted), 
             colour = subsidy) + 
      geom_point(size = 3,
                 aes(colour = subsidy),
                 position = position_dodge(0.5)) +
      geom_errorbar(aes(ymin = conf.low, 
                        ymax = conf.high,
                        colour = subsidy), 
                    width = 0.2,
                    position = position_dodge(0.5)) +
      geom_jitter(data = exp2_center,
                  mapping = aes(x = exposure,
                                y = get(parameter),
                                colour = subsidy),
                  alpha = 0.3) +
      scale_y_continuous(trans = "log10") +
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
            axis.line = element_line(colour = "black"))} 
  
  ## Another fork depending on which group of microorganism
  else if(type == "lines" & parameter == "bact_log_scale"){
    ## Make mini loop to frame depending on exposed or shaded microcosms
    for(i in 1:length(levels(exp2_center$exposure))){
    assign(paste0("ret", i),
           ggplot(data = model_effect %>%  
                     dplyr::filter(exposure == levels(exp2_center$exposure)[i]),
                   aes(x = get(parameter), 
                       y = predicted, 
                       colour = subsidy,
                       group = subsidy)) + 
            geom_line(aes(colour = subsidy)) +
            geom_ribbon(aes(ymin = conf.low, 
                            ymax = conf.high,
                            fill = subsidy), 
                        colour = NA,
                        alpha = 0.2) +
            geom_jitter(size = 3,
                       data = exp2_center %>% 
                         dplyr::filter(exposure == levels(exp2_center$exposure)[i]) ,
                       aes(colour = subsidy,
                           x = bact_log_scale,
                           y = chlorophyll_ugL),
                       height = 0.5) +
            ggtitle("") +
            xlab(xlab) +
            ylab(ylab) +
            xlim(-2, 2) +
            scale_y_continuous(trans = "log",
                               breaks = c(1, 20, 50, 400),
                               limits = c(0.1, 450)) +
            scale_color_manual(name = "Subsidy",
                               labels = c("Litter", "Litter + feces"), 
                               values = c("tan1", "tan4")) +
            scale_fill_manual(name = "Subsidy",
                               labels = c("Litter", "Litter + feces"), 
                               values = c("tan1", "tan4")) +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"))) }
    ## Join together
    ret <- 
      list(ret1, ret2)
  }
  else if(type == "lines" & parameter %in% c("chlorophyll_ugL_log_scale", "din_log_scale", "po4_log_scale")){
    ## Make mini loop to frame depending on exposed or shaded microcosms
    for(i in 1:length(levels(exp2_center$exposure))){
      assign(paste0("ret", i),
             ggplot(data = model_effect %>%  
                      dplyr::filter(exposure == levels(exp2_center$exposure)[i]),
                    aes(x = get(parameter), 
                        y = predicted, 
                        colour = subsidy,
                        group = subsidy)) + 
               geom_line(aes(colour = subsidy)) +
               geom_ribbon(aes(ymin = conf.low, 
                               ymax = conf.high,
                               fill = subsidy), 
                           colour = NA,
                           alpha = 0.2) +
               geom_jitter(size = 3,
                           data = exp2_center %>% 
                             dplyr::filter(exposure == levels(exp2_center$exposure)[i]) ,
                           aes(colour = subsidy,
                               x = get(parameter),
                               y = bact),
                           height = 0.5) +
               ggtitle("") +
               xlab(xlab) +
               ylab(ylab) +
               xlim(-2, 2) +
               scale_y_continuous(trans = "log",
                                   breaks = c(1, 2, 5, 10),
                                   limits = c(0.1, 10)) +
               scale_color_manual(name = "Subsidy",
                                  labels = c("Litter", "Litter + feces"), 
                                  values = c("tan1", "tan4")) +
               scale_fill_manual(name = "Subsidy",
                                 labels = c("Litter", "Litter + feces"), 
                                 values = c("tan1", "tan4")) +
               theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))) }
    ## Join together
    ret <- 
      list(ret1, ret2)
  }
    
    # Return
    return(ret)
  }
