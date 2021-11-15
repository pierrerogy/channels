# Load necessary packages
library(plater)
library(tidyverse)
library(readr)
library(stringr)
library(cowplot)
library(modelr)
library(lubridate)
library(broom)
source("R_code/functions.R")
library(assertr)


# Load necessary data
weeks <- 
  readr::read_csv("data_exp2/raw/weeks.csv") %>% 
  ## Unite date column for plate functions
  tidyr::unite(date, year, month, day, sep="", remove = FALSE) %>% 
  ## deal with the leading zero problem
  dplyr::mutate(date = ifelse(date == 2021111, 20211101, 
                       ifelse(date == 2021118, 20211108,
                              ifelse(date == 2021123, 20211203,
                                     date)))) %>% 
  dplyr::mutate(date = as.numeric(date))
chlorophyll <- 
  readr::read_csv("data_exp2/raw/chlorophyll.csv")
treatments <- 
  readr::read_csv("data_exp2/raw/treatments.csv")
mosquitoes <- 
  readr::read_csv("data_exp2/raw/mosquito_larvae.csv")
water <- 
  readr::read_csv("data_exp2/raw/water.csv")


# Most script generously provided by Kaleigh Davis (https://www.kaleighdavis.com/)
# Read in the data -- here we are making a list of datafiles
nutrient_files <-  
  c(list.files("data_exp2/raw/", 
               full.names = TRUE))
nut_parms <- 
  read.csv(file = "data_exp1/raw/nut_parms.csv")

# Sort files by nutrient
nox_files <- 
  nutrient_files[grepl("no", nutrient_files)]
phosphate_files <- 
  nutrient_files[grepl("po", nutrient_files)]
ammonium_files <- 
  nutrient_files[grepl("nh", nutrient_files)]
bacteria_files <- 
  nutrient_files[grepl("bac", nutrient_files)]


# Extract data from each file using custom function and remove empty wells
all_nox <- 
  nox_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))
all_pho <- 
  phosphate_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))
all_amm <- 
  ammonium_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))
all_bact <- 
  bacteria_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))

# Nitrite -----------------------------------------------------------------
# Prepare data
all_nitrite <- 
  all_nox %>%
  ## Make long form
  tidyr::gather(key = round, 
                value = absorbance, 
                no2_a, no2_b) %>% 
  ## Add plate column
  dplyr::mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  ## Add week
  dplyr::left_join(weeks %>% 
              dplyr::select(week, date)) %>% 
  ## Rename sample
  dplyr::rename(cup_number = no2no3_sam) %>% 
  ## Get rid of blank wells
  tidyr::unite(column, no2_stds, cup_number, 
               sep = "-", remove = T) %>% 
  dplyr::filter(column != "NA-NA") %>% 
  tidyr::separate(column, into = c("no2_stds", "cup_number"), 
                  sep = "-", remove = TRUE) %>% 
  dplyr::filter(!is.na(absorbance))

# Replace "NA" with actual NA
all_nitrite[all_nitrite == "NA"] = NA

# Check on structure of variables & make adjustments as necessary
str(all_nitrite)
all_nitrite <- 
  all_nitrite %>% 
  dplyr::mutate_at(vars(no2_stds), 
            ~as.numeric(.)) %>% 
  dplyr::mutate_at(vars(date, cup_number, week),
            ~as.character(.))

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
no2_split <- 
  all_nitrite %>%
  dplyr::group_by(week, plate) %>% 
  dplyr::group_split()

# Create data frame with zeroed cal curve values for each plate  
plate_zeroes_no2 <- 
  no2_split %>% 
  purrr::map_df(plate_specific_zeroes_no2)

# Small dataframe with only plate zero values
summ_plate_zeroes_no2 <- 
  plate_zeroes_no2 %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::summarise(mean_zero = mean(plate_zero))

# New column to sort by and nest dataset so it is compatible with the map function
nest_by_week_plate_no2 <- 
  plate_zeroes_no2 %>% 
  dplyr::group_by(week, plate) %>% 
  tidyr::nest() 

# Apply function to each plate to get data frame with all of your slopes and intercepts
nitrite_mod_parms <- 
  nest_by_week_plate_no2 %>% 
  dplyr::mutate(fit = map(data, no2_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(week, plate, tidied) %>% 
  tidyr::unnest(cols = c(tidied)) %>% 
  dplyr::select(week, plate, term, estimate) %>% 
  tidyr::spread(key = term, value = estimate) %>% 
  tidyr::gather(key = term, value = estimate, -week, -plate) %>% 
  dplyr::filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  dplyr::rename(intercept = `(Intercept)`,
         slope = no2_stds)

# Compare experimental values to expected ones
# Add plate zeroes to parameter table
nitrite_mod_parms1  <- 
  dplyr::left_join(nitrite_mod_parms, summ_plate_zeroes_no2)

# Getting ready to plop expected parameter values into this table
nit_parms <- 
  nut_parms %>% 
  dplyr::filter(nuts == "no2")

# Make vector for expected intercept to merge with modeled parameters table
nit_exp_int <- 
  rep(nit_parms$b, length = length(nitrite_mod_parms$slope))

# Make vector for expected slope to merge with modeled parameters table
nit_exp_slope <- 
  rep(nit_parms$m, length = length(nitrite_mod_parms$slope))

# Add expected values to the parm table
nitrite_mod_parms2 <- 
  cbind(nitrite_mod_parms1, nit_exp_int) %>% 
  dplyr::rename(nit_exp_int = ...6)
nitrite_mod_parms3 <- 
  cbind(nitrite_mod_parms2, nit_exp_slope) %>% 
  dplyr::rename(nit_exp_slope = ...7)

# Add columns for degree of divergence from expected vlaues
nitrite_mod_parms4 <- 
  nitrite_mod_parms3 %>% 
  dplyr::mutate(int_dist = abs(nit_exp_int - intercept),
         slope_dist = abs(nit_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
nitrite_mod_parms5 <- 
  nitrite_mod_parms4 %>% 
  dplyr::mutate(intercept = ifelse(int_dist < 0.4*nit_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept value
nitrite_mod_parms6 <- 
  nitrite_mod_parms5 %>% 
  dplyr::mutate(int_dist = nit_exp_int - intercept)

# If zero values are still too far from expected, throw out run
nit_dist_int <- 
  nitrite_mod_parms6 %>% 
  dplyr::filter(int_dist > 0.4*nit_exp_int)

# Look at slopes
nit_dist_slope <- 
  nitrite_mod_parms6 %>% 
  dplyr::filter(slope_dist > 0.4*nit_exp_slope)

# Predict concentrations with new parm values
nit_sams <- 
  all_nitrite %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::filter(!is.na(cup_number)) %>% 
  dplyr::select(-no2_stds)

# Calculate concentrations!
nit_sams1 <- 
  dplyr::left_join(nit_sams, nitrite_mod_parms6) %>% 
  dplyr::select(-nit_exp_int, -nit_exp_slope, -int_dist, -slope_dist) %>% 
  dplyr::mutate(no2_con = (absorbance - intercept) / slope) 

# Make negative concentrations and and positives below the LOD zeroes
no2_values <- 
  nit_sams1 %>% 
  dplyr::mutate(no2_con = ifelse(no2_con < 0, 0, no2_con)) %>%
  dplyr::select(-no2no3_a, -Wells, -no3_stds, -no2no3_b, 
                - absorbance, -intercept, -slope, -mean_zero) %>% 
  tidyr::pivot_wider(values_from = no2_con, names_from = plate) %>% 
  dplyr::rename(no2_a = A,
         no2_b = B)



# Nitrate -----------------------------------------------------------------
# Prepare data
all_nitrate <- 
  all_nox %>%
  ## Make long form
  tidyr::gather(key = round, 
         value = absorbance, 
         no2_a, no2_b, 
         no2no3_a, no2no3_b) %>% 
  ## Add plate column
  dplyr::mutate(replicate = ifelse(grepl("a$", round), 
                            "A", "B")) %>% 
  dplyr::mutate(Round = if_else(grepl("no3", round), 
                         "no2no3", "no2")) %>%
  dplyr::select(-round) %>% 
  dplyr::rename(round = Round) %>% 
  ## Get rid of blank wells
  dplyr::filter(!is.na(absorbance)) %>% 
  ## Spread data
  tidyr::spread(key = round, value = absorbance) %>% 
  # Add column for net absorbance due to nitrate
  dplyr::mutate(no3 = no2no3 - no2) %>% 
  dplyr::select(-Wells) %>% 
  ## Add week
  dplyr::left_join(weeks %>% 
                     dplyr::select(week, date)) %>% 
  ## Rename sample
  dplyr::rename(cup_number = no2no3_sam) 

# Some quick plotting to see how my calibration curves look
all_nitrate %>% 
  dplyr::filter(!is.na(no3_stds)) %>% 
  ggplot(aes(x = no3_stds, y = no3)) +
  geom_point() +
  facet_grid(week ~ replicate)

# Split data by plates
nat_split <- 
  all_nitrate %>%
  dplyr::group_by(week, replicate) %>% 
  dplyr::group_split()

# Make data frame with zeroed cal curve values for each plate  
full_plate_zeroes_no3 <- 
  nat_split %>% 
  purrr::map_df(plate_specific_zeroes_no3)

# Make small dataframe with only plate zero values
summ_plate_zeroes_no3 <- 
  full_plate_zeroes_no3 %>% 
  dplyr::group_by(week, replicate) %>% 
  dplyr::summarise(mean_zero = mean(plate_zero))

# New column to sort by and nest dataset so it is compatible with the map function
nest_by_week_rep_no3 <- 
  full_plate_zeroes_no3 %>% 
  dplyr::group_by(week, replicate) %>% 
  tidyr::nest() 

# Apply cal curve function to each plate to get data frame 
# with all slopes and intercepts
nitrate_mod_parms <- 
  nest_by_week_rep_no3 %>% 
  dplyr::mutate(fit = map(data, no3_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(week, replicate, tidied) %>% 
  tidyr::unnest(cols = c(tidied)) %>% 
  dplyr::select(week, replicate,  term, estimate) %>% 
  tidyr::spread(key = term, value = estimate) %>% 
  tidyr::gather(key = term, value = estimate, -week, -replicate) %>% 
  dplyr::filter(!is.na(estimate)) %>% 
  tidyr::spread(key= term, value = estimate) %>% 
  dplyr::rename(intercept = `(Intercept)`,
         slope = no3_stds)

# Compare experimental parm values to expected ones 
nitrate_mod_parms  <- 
  dplyr::left_join(nitrate_mod_parms, summ_plate_zeroes_no3)

# Getting ready to plop expected parameter values into this table
nat_parms <- 
  nut_parms %>% 
  dplyr::filter(nuts == "no3")

# Make vector for expected intercept to merge with modeled parameters table
nat_exp_int <- 
  rep(nat_parms$b, 
      length = length(nitrate_mod_parms$intercept))

# Make vector for expected slope to merge with modeled parameters table
nat_exp_slope <- 
  rep(nat_parms$m, 
      length = length(nitrate_mod_parms$intercept))

# Adding expected values to the parm table
nitrate_mod_parms <- 
  cbind(nitrate_mod_parms, nat_exp_int) %>% 
  dplyr::rename(nat_exp_int = ...6)
nitrate_mod_parms <- 
  cbind(nitrate_mod_parms, nat_exp_slope)%>% 
  dplyr::rename(nat_exp_slope = ...7)

# Add columns for degree of divergence from expected values
nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  dplyr::mutate(int_dist = abs(nat_exp_int - intercept),
         slope_dist = abs(nat_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  dplyr::mutate(intercept = ifelse(int_dist < 0.3*nat_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept values
nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  dplyr::mutate(int_dist = nat_exp_int - intercept)

# If zero values are still too far from expected, throw out run 
nat_dist_int <- 
  nitrate_mod_parms %>% 
  dplyr::filter(int_dist > 0.3*nat_exp_int)

# Look at slopes
nat_dist_slope <-
  nitrate_mod_parms %>% 
  dplyr::filter(abs(slope_dist) > 0.4*nat_exp_slope)

nat_sd_m <- 
  nitrate_mod_parms %>% 
  dplyr::group_by(week) %>% 
  dplyr::summarise(sd_m = sd(slope)) 

# Predict concentrations with new parm values 
nat_sams <- 
  all_nitrate %>% 
  dplyr::group_by(week, replicate) %>% 
  dplyr::filter(!is.na(cup_number)) %>% 
  dplyr::select(-no3_stds, -no2_stds)

#use lm parameters for each plate to estimate concentration
no3_values <- 
  dplyr::left_join(nat_sams, nitrate_mod_parms) %>% 
  dplyr::select(-nat_exp_int, -nat_exp_slope, -int_dist, -slope_dist) %>% 
  dplyr::mutate(no3_con = (no3 - intercept) / slope) %>% 
  dplyr::mutate(no3_con = ifelse(no3_con < 0, 0, no3_con)) %>% 
  dplyr::select(-no2no3, -no2, -no3, 
                -intercept, -slope, -mean_zero) %>% 
  tidyr::pivot_wider(values_from = no3_con, names_from = replicate) %>% 
  dplyr::rename(no3_a = A,
                no3_b = B)


# Phosphate ---------------------------------------------------------------
# Prepare data
all_phosphate <- 
  all_pho %>% 
  ## Make long form
  tidyr::gather(key = round, value = absorbance, po4_a, po4_b) %>%
  ## Add plate column
  dplyr::mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  ## Rename samples
  dplyr::rename(cup_number = po4_sam) %>% 
  ## Get rid of blank wells
  tidyr::unite(column, po4_stds, cup_number, sep = "-", remove = T) %>% 
  dplyr::filter(column != "NA-NA") %>% 
  tidyr::separate(column, into = c("po4_stds", "cup_number"), sep = "-", remove = TRUE) %>% 
  dplyr::filter(!is.na(absorbance)) %>%
  ## Add week
  dplyr::left_join(weeks %>% 
              dplyr::select(week, date))

#replace "NA" with actual NA
all_phosphate[all_phosphate == "NA"] = NA

# Check on structure of variables & make adjustments as necessary
str(all_phosphate)
all_phosphate <- 
  all_phosphate %>% 
  dplyr::mutate_at(vars(po4_stds), 
                   ~as.numeric(.)) %>% 
  dplyr::mutate_at(vars(date, cup_number, week), 
                   ~as.character(.))

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
po4_split <- 
  all_phosphate %>%
  dplyr::group_by(week, plate) %>% 
  dplyr::group_split()

# Create data frame with zeroed cal curve values for each plate  
plate_zeroes_po4 <- 
  po4_split %>% 
  purrr::map_df(plate_specific_zeroes_po4)

# Small dataframe with only plate zero values
summ_plate_zeroes_po4 <- 
  plate_zeroes_po4 %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::summarise(mean_zero = mean(plate_zero))

# New column to sort by and nest dataset so it is compatible with the map function
nest_by_week_plate_po4 <- 
  plate_zeroes_po4 %>% 
  dplyr::group_by(week, plate) %>% 
  tidyr::nest() 

# Apply function to each plate to get data frame with all of your slopes and intercepts
phosphate_mod_parms <- 
  nest_by_week_plate_po4 %>% 
  dplyr::mutate(fit = map(data, po4_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(week, plate, tidied) %>% 
  tidyr::unnest(cols = c(tidied)) %>% 
  dplyr::select(week, plate, term, estimate) %>% 
  tidyr::spread(key = term, value = estimate) %>% 
  tidyr::gather(key = term, value = estimate, -week, -plate) %>% 
  dplyr::filter(!is.na(estimate)) %>% 
  tidyr::spread(key= term, value = estimate) %>% 
  dplyr::rename(intercept = `(Intercept)`,
         slope = po4_stds)

# Compare experimental values to expected ones
# Add plate zeroes to parameter table
phosphate_mod_parms1  <- 
  dplyr::left_join(phosphate_mod_parms, summ_plate_zeroes_po4)

# Getting ready to plop expected parameter values into this table
pho_parms <- 
  nut_parms %>% 
  dplyr::filter(nuts == "po4")

# Make vector for expected intercept to merge with modeled parameters table
pho_exp_int <- 
  rep(pho_parms$b, length = length(phosphate_mod_parms$slope))

# Make vector for expected slope to merge with modeled parameters table
pho_exp_slope <- 
  rep(pho_parms$m, length = length(phosphate_mod_parms$slope))

# Add expected values to the parm table
phosphate_mod_parms2 <- 
  cbind(phosphate_mod_parms1, pho_exp_int) %>% 
  dplyr::rename(pho_exp_int = ...6)
phosphate_mod_parms3 <- 
  cbind(phosphate_mod_parms2, pho_exp_slope) %>% 
  dplyr::rename(pho_exp_slope = ...7)

# Add columns for degree of divergence from expected vlaues
phosphate_mod_parms4 <- 
  phosphate_mod_parms3 %>% 
  dplyr::mutate(int_dist = abs(pho_exp_int - intercept),
         slope_dist = abs(pho_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
phosphate_mod_parms5 <- 
  phosphate_mod_parms4 %>% 
  dplyr::mutate(intercept = ifelse(int_dist < 0.4*pho_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept value
phosphate_mod_parms6 <- 
  phosphate_mod_parms5 %>% 
  dplyr::mutate(int_dist = pho_exp_int - intercept)

# If zero values are still too far from expected, throw out run
pho_dist_int <- 
  phosphate_mod_parms6 %>% 
  dplyr::filter(int_dist > 0.4*pho_exp_int)

# Look at slopes
pho_dist_slope <- 
  phosphate_mod_parms6 %>% 
  dplyr::filter(slope_dist > 0.4*pho_exp_slope)

# Predict concentrations with new parm values
pho_sams <- 
  all_phosphate %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::filter(!is.na(cup_number)) %>% 
  dplyr::select(-po4_stds)

# Calculate concentrations!
pho_sams1 <- 
  dplyr::left_join(pho_sams, phosphate_mod_parms6) %>% 
  dplyr::select(-pho_exp_int, -pho_exp_slope, -int_dist, -slope_dist) %>% 
  dplyr:: mutate(po4_con = (absorbance - intercept) / slope) 

# Make negative concentrations zeroes
po4_values <- 
  pho_sams1 %>% 
  dplyr::mutate(po4_con = ifelse(po4_con < 0, 0, po4_con)) %>% 
  dplyr::select(-slope, -mean_zero, -intercept, -absorbance, -Wells) %>%
  tidyr::pivot_wider(values_from = po4_con, names_from = plate) %>% 
  dplyr::rename(po4_a = A,
         po4_b = B)


# Ammonium ----------------------------------------------------------------
# Prepare data
all_ammonium <- 
  all_amm %>% 
  ## Make long form
  tidyr::gather(key = round, value = absorbance, nh4_a, nh4_b) %>%
  ## Add plate column
  dplyr::mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  ## Rename samples
  dplyr::rename(cup_number = nh4_sam) %>% 
  ## Get rid of blank wells
  tidyr::unite(column, nh4_stds, cup_number, sep = "-", remove = T) %>% 
  dplyr::filter(column != "NA-NA") %>% 
  tidyr::separate(column, into = c("nh4_stds", "cup_number"), sep = "-", remove = TRUE) %>% 
  dplyr::filter(!is.na(absorbance)) %>%
  ## Add week
  dplyr::left_join(weeks %>% 
                     dplyr::select(week, date))

# Replace "NA" with actual NA
all_ammonium[all_ammonium =="NA"] = NA

# C heck on structure of variables & make adjustments as necessary
str(all_ammonium)
all_ammonium <- 
  all_ammonium %>% 
  dplyr::mutate_at(vars(nh4_stds), ~as.numeric(.)) %>% 
  dplyr::mutate_at(vars(date, cup_number, week), ~as.character(.)) 

# Some quick plotting to see how my calibration curves look
all_ammonium %>% 
  dplyr::filter(!is.na(nh4_stds)) %>% 
  ggplot(aes(x = nh4_stds, y = absorbance)) +
  geom_point() +
  facet_grid(week ~ plate)

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
nh4_split <- 
  all_ammonium %>%
  dplyr::group_by(week, plate) %>% 
  dplyr::group_split()

# Make data frame with zeroed cal curve values for each plate  
plate_zeroes_nh4 <- 
  nh4_split %>% 
  purrr::map_df(plate_specific_zeroes_nh4)

# Get plate values
summ_plate_zeroes_nh4 <- 
  plate_zeroes_nh4 %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::summarise(mean_zero = mean(plate_zero))


# New column to sort by and nest dataset so it is compatible with the map function
nest_by_date_plate_nh4 <- 
  plate_zeroes_nh4 %>% 
  dplyr::group_by(week, plate) %>% 
  tidyr::nest() 


# Apply function to each plate to get data frame with each slope and intercept
ammonium_mod_parms <- 
  nest_by_date_plate_nh4 %>% 
  dplyr::mutate(fit = map(data, nh4_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(week, plate, tidied) %>% 
  tidyr::unnest(cols = c(tidied)) %>% 
  dplyr::select(week, plate, term, estimate) %>% 
  tidyr::spread(key = term, value = estimate) %>% 
  tidyr::gather(key = term, value = estimate, -week, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  tidyr::spread(key= term, value = estimate) %>% 
  dplyr::rename(intercept = `(Intercept)`,
         slope = nh4_stds)

# Compare experimental parm values to expected ones 
# Add plate zeroes to parameter table
ammonium_mod_parms  <- 
  dplyr::left_join(ammonium_mod_parms, summ_plate_zeroes_nh4)

# Getting ready to plop expected parameter values into this table
amm_parms <- 
  nut_parms %>% 
  dplyr::filter(nuts == "nh4")

# Make vector for expected intercept to merge with modeled parameters table
amm_exp_int <- 
  rep(amm_parms$b, length = length(ammonium_mod_parms$slope))

# Make vector for expected slope to merge with modeled parameters table
amm_exp_slope <- 
  rep(amm_parms$m, length = length(ammonium_mod_parms$slope))

# Add expected values to the parm table
ammonium_mod_parms <- 
  cbind(ammonium_mod_parms, amm_exp_int) %>% 
  rename(amm_exp_int = ...6)
ammonium_mod_parms <- 
  cbind(ammonium_mod_parms, amm_exp_slope) %>% 
  rename(amm_exp_slope = ...7)

# Add columns for degree of divergence from expected vlaues
ammonium_mod_parms <- 
  ammonium_mod_parms %>% 
  mutate(int_dist = abs(amm_exp_int - intercept),
         slope_dist = abs(amm_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
ammonium_mod_parms <- 
  ammonium_mod_parms %>% 
  mutate(intercept = ifelse(int_dist < 0.4*amm_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept value
ammonium_mod_parms <- 
  ammonium_mod_parms %>% 
  mutate(int_dist = amm_exp_int - intercept)

# If zero values are still too far from expected, throw out run
nh4_dist_int <- 
  ammonium_mod_parms %>% 
  filter(int_dist > 0.4*amm_exp_int)

# Look at slopes
nh4_dist_slope <- 
  ammonium_mod_parms %>% 
  filter(slope_dist > 0.4*amm_exp_slope) 

# Predict concentrations with new parm values 
amm_split <- 
  all_ammonium %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::filter(!is.na(cup_number)) %>% 
  dplyr::select(-nh4_stds)

# Calculate concentrations!
amm_split <- 
  dplyr::left_join(amm_split, ammonium_mod_parms) %>% 
  dplyr::select(-amm_exp_int, -amm_exp_slope, -int_dist, -slope_dist, -Wells) %>% 
  dplyr::mutate(nh4_con = (absorbance - intercept) / slope)

# Negative concentrations = 0
nh4_values <- 
  amm_split %>% 
  mutate(nh4_con = ifelse(nh4_con < 0, 0, nh4_con)) %>% 
  dplyr::select(-mean_zero, -absorbance, -intercept, -slope) %>% 
  tidyr::pivot_wider(values_from = nh4_con, names_from = plate) %>% 
  dplyr::rename(nh4_a = A,
         nh4_b = B)



# Bacteria ----------------------------------------------------------------
# Add in replicate information and actual concentration of spheres in standard
# Prepare data
all_bacteria <- 
  all_bact %>% 
  ## Make long form
  tidyr::gather(key = round, value = absorbance, bact_a, bact_b) %>%
  ## Add plate column
  dplyr::mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  ## Rename samples
  dplyr::rename(cup_number = bact_sam) %>% 
  ## Get rid of blank wells
  tidyr::unite(column, bact_stds, cup_number, sep = "-", remove = T) %>% 
  dplyr::filter(column != "NA-NA") %>% 
  tidyr::separate(column, into = c("bact_stds", "cup_number"), sep = "-", remove = TRUE) %>% 
  dplyr::filter(!is.na(absorbance)) %>%
  ## Add week
  dplyr::left_join(weeks %>% 
                     dplyr::select(week, date))

# Replace "NA" with actual NA
all_bacteria[all_bacteria =="NA"] = NA

# C heck on structure of variables & make adjustments as necessary
str(all_bacteria)
all_bacteria <- 
  all_bacteria %>% 
  dplyr::mutate_at(vars(bact_stds), 
                   ~as.numeric(.)) %>% 
  dplyr:: mutate_at(vars(date, cup_number, week), 
                    ~as.character(.)) 

# Some quick plotting to see how my calibration curves look
all_bacteria %>% 
  filter(!is.na(bact_stds)) %>% 
  ggplot(aes(x = bact_stds, y = absorbance)) +
  geom_point() +
  facet_grid(week ~ plate) 

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
bact_split <- 
  all_bacteria %>%
  dplyr::group_by(week, plate) %>% 
  dplyr::group_split()

# Make data frame with zeroed cal curve values for each plate  
plate_zeroes_bact <- 
  bact_split %>% 
  purrr::map_df(plate_specific_zeroes_bact)

# Get plate values
summ_plate_zeroes_bact <- 
  plate_zeroes_bact %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::summarise(mean_zero = mean(plate_zero))


# New column to sort by and nest dataset so it is compatible with the map function
nest_by_date_plate_bact <- 
  plate_zeroes_bact %>% 
  dplyr::group_by(week, plate) %>% 
  tidyr::nest() 


# Apply function to each plate to get data frame with each slope and intercept
bacteria_mod_parms <- 
  nest_by_date_plate_bact %>% 
  dplyr::mutate(fit = map(data, bact_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(week, plate, tidied) %>% 
  tidyr::unnest(cols = c(tidied)) %>% 
  dplyr::select(week, plate, term, estimate) %>% 
  tidyr::spread(key = term, value = estimate) %>% 
  tidyr::gather(key = term, value = estimate, -week, -plate) %>% 
  dplyr::filter(!is.na(estimate)) %>% 
  tidyr::spread(key= term, value = estimate) %>% 
  dplyr::rename(intercept = `(Intercept)`,
         slope = bact_stds)


# Predict concentrations with new parm values 
bact_split <- 
  all_bacteria %>% 
  dplyr::group_by(week, plate) %>% 
  dplyr::filter(!is.na(cup_number)) %>% 
  dplyr::select(-bact_stds)

# Calculate concentrations!
bact_split <- 
  left_join(bact_split, bacteria_mod_parms) %>% 
  mutate(bact_con = (absorbance - intercept) / slope)



# Negative concentrations = 0
## First visit NA because different, uncomparable method
bact_values <- 
  bact_split %>% 
  dplyr::mutate(bact_con = ifelse(bact_con < 0,
                                  0, bact_con)) %>%
  dplyr::select(-Wells, - absorbance, -intercept, -slope) %>% 
  tidyr::pivot_wider(values_from = bact_con, names_from = plate) %>% 
  dplyr::rename(bact_a = A,
                bact_b = B)


# Convert chlorophyll values ----------------------------------------------
chlorophyll <- 
  chlorophyll %>% 
  ## Remove blank values
  dplyr::left_join(chlorophyll %>% 
              filter(cup_number == "acetone") %>% 
              dplyr::select(week, chlorophyll_raw) %>% 
              group_by(week) %>% 
              summarise_all(mean) %>% 
              rename(blanks = chlorophyll_raw)) %>% 
  dplyr::mutate(chlorophyll_raw = chlorophyll_raw - blanks) %>% 
  dplyr::filter(cup_number != "acetone") %>% 
  ## Get concentration!
  dplyr::mutate(chlorophyll_ugL = chlorophyll_raw * (solvent_added / v_filtered)) %>% 
  dplyr::select(-v_filtered, -solvent_added, -v_solvent_entered,
                -v_water_entered, -chlorophyll_raw, -blanks)



# Combine and save data ---------------------------------------------------
weekly_data <- 
  water %>%
  dplyr::select(-day) %>% 
  dplyr::mutate(week = as.character(week),
         cup_number = as.character(cup_number)) %>% 
  dplyr::left_join(treatments %>% 
                     dplyr::mutate(cup_number = as.character(cup_number))) %>% 
  dplyr::left_join(chlorophyll %>% 
              dplyr::mutate(week = as.character(week),
                     cup_number = as.character(cup_number)),
              by = c("week", "cup_number")) %>% 
  left_join(po4_values,
            by = c("week", "cup_number")) %>% 
  left_join(nh4_values %>% 
              dplyr::select(-date),
            by = c("week", "cup_number")) %>% 
  left_join(no2_values %>% 
              dplyr::select(-date),
            by = c("week", "cup_number")) %>% 
  left_join(no3_values %>% 
              dplyr::select(-date),
            by = c("week", "cup_number")) %>% 
  left_join(bact_values %>% 
              dplyr::select(-date),
            by = c("week", "cup_number")) %>% 
  dplyr::mutate(date = lubridate::as_date(date)) %>% 
  ## Average two values
  dplyr::rowwise() %>% 
  dplyr::mutate(po4 = mean(c(po4_a, po4_b), na.rm = T),
                no2 = mean(c(no2_a, no2_b), na.rm = T),
                nh4 = mean(c(nh4_a, nh4_b), na.rm = T),
                no3 = mean(c(no3_a, no3_b), na.rm = T),
                bact = mean(c(bact_a, bact_b), na.rm = T)) %>% 
  dplyr::select(-po4_a, -po4_b, -nh4_a, -nh4_b,
                -no2_a, -no2_b, -no3_a, -no3_b,
                -bact_a, -bact_b)

# Save data
readr::write_csv(weekly_data,
                 "data_exp2/weekly_measurements_exp2.csv")
readr::write_csv(weekly_data,
                 "microchannels/appdata/weekly_measurements_exp2.csv")


# Prepare mosquito data ---------------------------------------------------------------------
# Clean data
mosquitoes_tidied <- 
  mosquitoes %>% 
  dplyr::mutate(date = lubridate::as_date(lubridate::dmy(Day))) %>% 
  ## Add treatments
  dplyr::left_join(treatments,
                   by = "cup_number") %>%
  ## make missing dead
  dplyr::mutate(event = ifelse(event == "missing",
                               "death", event))
  

# Save data
readr::write_csv(mosquitoes_tidied,
                 "microchannels/appdata/mosquitoes.csv")
readr::write_csv(mosquitoes_tidied,
                 "data_exp2/mosquitoes.csv")

# Assertions --------------------------------------------------------------


# Check data format
str(weekly_measurements)

# Make sure all data are in OK range
weekly_measurements %>% 
  mutate(n_leaves = as.numeric(n_leaves),
         larvae = as.numeric(larvae),
         subsidy_id = as.numeric(subsidy_id),
         capacity_mL = as.numeric(capacity_mL)) %>% 
  dplyr::select(-vessel, -shading, -cup_number, -week) %>% 
  assert(within_bounds(0, Inf),
         everything()) %>% 
  verify(nrow(.) == 312)

# Make sure no typos in several columns
## 12 visits
length(unique(weekly_measurements$week)) == 12
## 52 bromelias
length(unique(weekly_measurements$cup_number)) == 52
## 52 subsidies
length(unique(weekly_measurements$subsidy_id)) == 52
## Two x two level treatments
length(unique(weekly_measurements$vessel)) == 2
length(unique(weekly_measurements$shading)) == 2


# Make a bunch of histograms
## First select series of columns to plot
index_vec <- 
  c(3, 4, 6, 10:27)
## Make a loop for all plots
par(mfrow = c(2,2))
for(i in index_vec){
  hist(weekly_measurements[,i],
       main = paste0(colnames(weekly_measurements[i])))
}
par(mfrow = c(1,1))

