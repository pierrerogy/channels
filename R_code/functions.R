# Load necessary packages
# Functions written by Kaleigh Davis (https://www.kaleighdavis.com/)


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
# Get coefficients from model output
get_coefs <- function(data, model, bacteria = FALSE, chloro = FALSE, feces = FALSE){
  
  ## Extract coefficients if bacteria and no algae
  if(bacteria == TRUE){
    coefs <- 
      data %>% 
      dplyr::select(shading, din_scale, bact_scale) %>% 
      dplyr::bind_cols(predict(model, 
                               marginal = model$Random$formula,
                               type = "response",
                               interval = "confidence") %>% 
                         tibble::as_tibble() %>% 
                         dplyr::mutate_all(exp))
  }
  ## Extract coefficients if no bacteria but algae
  if(chloro == TRUE){
    coefs <- 
      data %>% 
      dplyr::select(shading, din_scale, chloro_scale) %>% 
      dplyr::bind_cols(predict(model, 
                               marginal = model$Random$formula,
                               type = "response",
                               interval = "confidence") %>% 
                         tibble::as_tibble() %>% 
                         dplyr::mutate_all(exp))
      
  }
  ## Extract coefficients just for feces
  if(feces == TRUE){
    ## Get coefficients
    coefs <- 
      data %>% 
      dplyr::select(shading, subsidy) %>% 
      dplyr::bind_cols(predict(model, 
                               marginal = model$Random$formula,
                               type = "response",
                               interval = "confidence") %>% 
                         tibble::as_tibble() %>% 
                         dplyr::mutate_all(exp)) %>% 
      dplyr::group_by(shading, subsidy) %>% 
      dplyr::summarise_all(mean)
  }

  ## Return
  return(coefs)
}

# Build point coordinates for given predictor, and add confint
get_points <- function(data, coefs){
  ## Empty data frame to return
  ret <- 
    data.frame()
  
  ## Loop going through parameters
  for(parameter in unique(coefs$what)){
  ## Get minimum and maximum values for parameter
  min <- 
    data %>% 
    dplyr::select(paste0(parameter, "_scale")) %>% 
    min(na.rm = TRUE)
  max <- 
    data %>% 
    dplyr::select(paste0(parameter, "_scale")) %>% 
    max(na.rm = TRUE)
  
  ## Make vector of points between this min and max, every 0.2
  values <- 
    data.frame(predictor = seq(from = min,
                               to = max,
                               by = 0.2))
  
  ## Get intercept and slope from coefficients
  equation <- 
    coefs %>% 
    dplyr::filter(what == parameter)
  
  ## Get points for each intercept and slope
  for(i in 1:nrow(equation)){
    ## Get equation parameters
    row <- 
      equation[i,]
    intercept <- 
      row[3] %>% 
      pull()
    slope <- 
      row[4] %>% 
      pull()
    ci <- 
      row[5] %>% 
      pull()
    ## Get prediction for each point of predictor
    vals <- 
      values %>% 
      dplyr::mutate(predicted = intercept + predictor*slope,
                    ci_low = predicted - ci,
                    ci_high = predicted + ci,
                    what = pull(row[1]),
                    shading = pull(row[2]))
    
    ## Stick to ret
    ret <- 
      ret %>% 
      dplyr::bind_rows(vals)
    }
    
   
      
    
    
  }
  
  ## Return
  return(ret)
  
}
  
