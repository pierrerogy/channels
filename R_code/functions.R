# Load necessary packages
# Functions written by Kaleigh Davis (https://www.kaleighdavis.com/)

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
