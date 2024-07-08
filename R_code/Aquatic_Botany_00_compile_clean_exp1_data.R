# Load necessary packages
library(plater)
library(tidyverse)
library(stringr)
library(cowplot)
library(modelr)
library(lubridate)
library(broom)
source("R_code/functions.R")
library(assertr)


# Load necessary data
visits <- 
  read.csv("data_exp1/raw/visits.csv") %>% 
  ## Unite date column for plate functions
  unite(date, year, month, day, sep="", remove = FALSE) %>% 
  ## deal with the leading zero problem
  mutate(date = ifelse(date == 2020115, 20201105, 
                       ifelse(date == 2020119, 20201109,
                              ifelse(date == 2020128, 20201208,
                                     ifelse(date == 2020129, 20201209,
                              date))))) %>% 
  mutate(date = as.numeric(date))
weekly_measurements_init <- 
  read.csv("data_exp1/raw/weekly_measurements.csv")
chlorophyll <- 
  read.csv("data_exp1/raw/chlorophyll.csv")[,1:8]
bromeliads <- 
  read.csv("data_exp1/raw/bromeliads.csv")
subsidies <- 
  read.csv("data_exp1/raw/subsidies.csv")


# Most script generously provided by Kaleigh Davis (https://www.kaleighdavis.com/)
# Read in the data -- here we are making a list of datafiles
nutrient_files <-  
  c(list.files("data_exp1/raw/", 
               full.names = TRUE))
nut_parms <- 
  read.csv(file = "data_exp1/raw/nut_parms.csv")

# Sort files by nutrient
nox_files <- 
  nutrient_files[grepl("_no", nutrient_files)]
phosphate_files <- 
  nutrient_files[grepl("_po", nutrient_files)]
ammonium_files <- 
  nutrient_files[grepl("_nh", nutrient_files)]
bacteria_files <- 
  nutrient_files[grepl("_bac", nutrient_files)]


# Extract data from each file using custom function and remove empty wells
all_nox <- 
  nox_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))
all_phosphate <- 
  phosphate_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))
all_ammonium <- 
  ammonium_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))
all_bacteria <- 
  bacteria_files %>% 
  map_df(read_plate_function) %>% 
  filter(!is.na(date))

# Nitrite -----------------------------------------------------------------
# Make data long form, add plate column
all_nitrite2 <- 
  all_nox %>% 
  gather(key = round, value = absorbance, no2_a, no2_b) %>% 
  mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  left_join(visits %>% 
              dplyr::select(visit_id, date)) %>% 
  mutate(visit_id = ifelse(grepl("H0_", no2no3_sam),
                           "H0", 
                           ifelse(grepl("H1_", no2no3_sam),
                                  "H1",
                                  ifelse(grepl("H4_", no2no3_sam),
                                         "H4",
                                         ifelse(grepl("H8_", no2no3_sam),
                                                "H8",
                                                visit_id))))) %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H0_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H1_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H4_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H8_",
            replacement = "")

# Isolate rows with combined visit id
new_rows <- 
  all_nitrite2 %>% 
  filter(visit_id == "H0H1H4H8")

# Remove rows from dataframe and add mutated ones
all_nitrite2 <-
  all_nitrite2 %>% 
  anti_join(new_rows) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H0")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H1")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H4")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H8"))

# Get rid of blank wells
all_nitrite2 <- 
  all_nitrite2 %>% 
  rename(bromeliad_id = no2no3_sam) %>% 
  unite(column, no2_stds, bromeliad_id, sep = "-", remove = T) %>% 
  filter(column != "NA-NA") %>% 
  separate(column, into = c("no2_stds", "bromeliad_id"), sep = "-", remove = TRUE)

# Replace "NA" with actual NA
all_nitrite2[all_nitrite2 == "NA"] = NA

# Check on structure of variables & make adjustments as necessary
str(all_nitrite2)
all_nitrite2 <- 
  all_nitrite2 %>% 
  mutate_at(vars(no2_stds), ~as.numeric(.)) %>% 
  mutate_at(vars(date, bromeliad_id, visit_id), ~as.character(.))

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
no2_split <- 
  all_nitrite2 %>%
  group_by(visit_id, plate) %>% 
  group_split()

# Create data frame with zeroed cal curve values for each plate  
plate_zeroes_no2 <- 
  no2_split %>% 
  map_df(plate_specific_zeroes_no2)

# Small dataframe with only plate zero values
summ_plate_zeroes_no2 <- 
  plate_zeroes_no2 %>% 
  group_by(visit_id, plate) %>% 
  summarise(mean_zero = mean(plate_zero))

# New column to sort by and nest dataset so it is compatible with the map function
nest_by_week_plate_no2 <- 
  plate_zeroes_no2 %>% 
  group_by(visit_id, plate) %>% 
  nest() 

# Apply function to each plate to get data frame with all of your slopes and intercepts
nitrite_mod_parms <- 
  nest_by_week_plate_no2 %>% 
  mutate(fit = map(data, no2_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(visit_id, plate, tidied) %>% 
  unnest(cols = c(tidied)) %>% 
  dplyr::select(visit_id, plate, term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -visit_id, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = no2_stds)

# Compare experimental values to expected ones
# Add plate zeroes to parameter table
nitrite_mod_parms1  <- 
  left_join(nitrite_mod_parms, summ_plate_zeroes_no2)

# Getting ready to plop expected parameter values into this table
nit_parms <- 
  nut_parms %>% 
  filter(nuts == "no2")

# Make vector for expected intercept to merge with modeled parameters table
nit_exp_int <- 
  rep(nit_parms$b, length = length(nitrite_mod_parms$slope))

# Make vector for expected slope to merge with modeled parameters table
nit_exp_slope <- 
  rep(nit_parms$m, length = length(nitrite_mod_parms$slope))

# Add expected values to the parm table
nitrite_mod_parms2 <- 
  cbind(nitrite_mod_parms1, nit_exp_int) %>% 
  rename(nit_exp_int = ...6)
nitrite_mod_parms3 <- 
  cbind(nitrite_mod_parms2, nit_exp_slope) %>% 
  rename(nit_exp_slope = ...7)

# Add columns for degree of divergence from expected vlaues
nitrite_mod_parms4 <- 
  nitrite_mod_parms3 %>% 
  mutate(int_dist = abs(nit_exp_int - intercept),
         slope_dist = abs(nit_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
nitrite_mod_parms5 <- 
  nitrite_mod_parms4 %>% 
  mutate(intercept = ifelse(int_dist < 0.4*nit_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept value
nitrite_mod_parms6 <- 
  nitrite_mod_parms5 %>% 
  mutate(int_dist = nit_exp_int - intercept)

# If zero values are still too far from expected, throw out run
nit_dist_int <- 
  nitrite_mod_parms6 %>% 
  filter(int_dist > 0.4*nit_exp_int)

# Look at slopes
nit_dist_slope <- 
  nitrite_mod_parms6 %>% 
  filter(slope_dist > 0.4*nit_exp_slope)

# Predict concentrations with new parm values
nit_sams <- 
  all_nitrite2 %>% 
  group_by(visit_id, plate) %>% 
  filter(!is.na(bromeliad_id)) %>% 
  dplyr::select(-no2_stds)

# Calculate concentrations!
nit_sams1 <- 
  left_join(nit_sams, nitrite_mod_parms6) %>% 
  dplyr::select(-nit_exp_int, -nit_exp_slope, -int_dist, -slope_dist) %>% 
  mutate(no2_con = (absorbance - intercept) / slope) 

# Calculate LOD
#LOD = 3*sd(slope)/slope
nit_sd_m <- 
  nitrite_mod_parms6 %>% 
  group_by(visit_id) %>% 
  summarise(sd_m = sd(slope))

nit_mean_m <- 
  nitrite_mod_parms6 %>% 
  group_by(visit_id) %>%
  summarise(mean_m = mean(slope)) %>% 
  left_join(nit_sd_m) %>% 
  mutate(no2_lod = 3*(sd_m)/(mean_m))

# Make negative concentrations and and positives below the LOD zeroes
no2_values <- 
  nit_sams1 %>% 
  left_join(nit_mean_m %>% 
              dplyr::select(visit_id, no2_lod)) %>% 
  mutate(no2_con = ifelse(no2_con < 0 | no2_con < no2_lod, 0, no2_con)) %>% 
  pivot_wider(values_from = no2_con, names_from = plate) %>% 
  rename(no2_a = A,
         no2_b = B)
  


# Nitrate -----------------------------------------------------------------
# Make data long form, add plate column, change round to just be no2 or no2no3
all_nitrate2 <- 
  all_nox %>% 
  gather(key = round, 
         value = absorbance, 
         no2_a, no2_b, no2no3_a, no2no3_b) %>% 
  mutate(replicate = ifelse(grepl("a$", round), 
                            "A", "B")) %>% 
  mutate(Round = if_else(grepl("no3", round), 
                         "no2no3", "no2")) %>%
  dplyr::select(-round) %>% 
  rename(round = Round) %>% 
  filter(!is.na(absorbance)) %>% 
  spread(key = round, value = absorbance) %>% 
  #add column for net absorbance due to nitrate
  mutate(no3 = no2no3 - no2) %>% 
  left_join(visits %>% 
              dplyr::select(date, visit_id)) %>% 
  dplyr::select(-Wells) %>% 
  rename(bromeliad_id = no2no3_sam) %>% 
  mutate(visit_id = ifelse(grepl("H0_", bromeliad_id),
                           "H0", 
                           ifelse(grepl("H1_", bromeliad_id),
                                  "H1",
                                  ifelse(grepl("H4_", bromeliad_id),
                                         "H4",
                                         ifelse(grepl("H8_", bromeliad_id),
                                                "H8",
                                                visit_id))))) %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H0_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H1_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H4_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H8_",
            replacement = "")

# Isolate
new_rows <- 
  all_nitrate2 %>% 
  filter(visit_id == "H0H1H4H8")

# Remove rows from dataframe and add mutated ones
all_nitrate2 <-
  all_nitrate2 %>% 
  anti_join(new_rows) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H0")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H1")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H4")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H8"))

# Some quick plotting to see how my calibration curves look
all_nitrate2 %>% 
  filter(!is.na(no3_stds)) %>% 
  ggplot(aes(x = no3_stds, y = no3)) +
  geom_point() +
  facet_grid(visit_id ~ replicate)

# Split data by plates
nat_split <- 
  all_nitrate2 %>%
  group_by(visit_id, replicate) %>% 
  group_split()

# Make data frame with zeroed cal curve values for each plate  
full_plate_zeroes_no3 <- 
  nat_split %>% 
  map_df(plate_specific_zeroes_no3)

# Make small dataframe with only plate zero values
summ_plate_zeroes_no3 <- 
  full_plate_zeroes_no3 %>% 
  group_by(visit_id, replicate) %>% 
  summarise(mean_zero = mean(plate_zero))

# New column to sort by and nest dataset so it is compatible with the map function
nest_by_week_rep_no3 <- 
  full_plate_zeroes_no3 %>% 
  group_by(visit_id, replicate) %>% 
  nest() 

# Apply cal curve function to each plate to get data frame 
# with all slopes and intercepts
nitrate_mod_parms <- 
  nest_by_week_rep_no3 %>% 
  mutate(fit = map(data, no3_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(visit_id, replicate, tidied) %>% 
  unnest(cols = c(tidied)) %>% 
  dplyr::select(visit_id, replicate,  term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -visit_id, -replicate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = no3_stds)

# Compare experimental parm values to expected ones 
nitrate_mod_parms  <- 
  left_join(nitrate_mod_parms, summ_plate_zeroes_no3)

# Getting ready to plop expected parameter values into this table
nat_parms <- 
  nut_parms %>% 
  filter(nuts == "no3")

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
  rename(nat_exp_int = ...6)
nitrate_mod_parms <- 
  cbind(nitrate_mod_parms, nat_exp_slope)%>% 
  rename(nat_exp_slope = ...7)

# Add columns for degree of divergence from expected values
nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  mutate(int_dist = abs(nat_exp_int - intercept),
         slope_dist = abs(nat_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  mutate(intercept = ifelse(int_dist < 0.3*nat_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept values
nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  mutate(int_dist = nat_exp_int - intercept)

# If zero values are still too far from expected, throw out run 
nat_dist_int <- 
  nitrate_mod_parms %>% 
  filter(int_dist > 0.3*nat_exp_int)

# Look at slopes
nat_dist_slope <-
  nitrate_mod_parms %>% 
  filter(abs(slope_dist) > 0.4*nat_exp_slope)

nat_sd_m <- 
nitrate_mod_parms %>% 
  group_by(visit_id) %>% 
  summarise(sd_m = sd(slope)) 

nat_mean_m <- 
  nitrate_mod_parms %>% 
  group_by(visit_id) %>%
  summarise(mean_m = mean(slope)) %>% 
  left_join(nat_sd_m) %>% 
  mutate(no3_lod = 3*(sd_m)/(mean_m))

nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  left_join(nat_mean_m %>% 
              dplyr::select(visit_id, no3_lod))

# Predict concentrations with new parm values 
nat_sams <- 
  all_nitrate2 %>% 
  group_by(visit_id, replicate) %>% 
  filter(!is.na(bromeliad_id)) %>% 
  dplyr::select(-no3_stds, -no2_stds)

#use lm parameters for each plate to estimate concentration
no3_values <- 
  left_join(nat_sams, nitrate_mod_parms) %>% 
  dplyr::select(-nat_exp_int, -nat_exp_slope, -int_dist, -slope_dist) %>% 
  mutate(no3_con = (no3 - intercept) / slope) %>% 
  mutate(no3_con = ifelse(no3_con < 0 | no3_con < no3_lod, 
                          0, no3_con)) %>% 
  pivot_wider(values_from = no3_con, names_from = replicate) %>% 
  rename(no3_a = A,
         no3_b = B)


# Phosphate ---------------------------------------------------------------
# Make data long form, add plate column
all_phosphate2 <- 
  all_phosphate %>% 
  gather(key = round, value = absorbance, po4_a, po4_b) %>% 
  mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  left_join(visits %>% 
              dplyr::select(visit_id, date))%>% 
  mutate(visit_id = ifelse(grepl("H0_", po4_sam),
                           "H0", 
                           ifelse(grepl("H1_", po4_sam),
                                  "H1",
                                  ifelse(grepl("H4_", po4_sam),
                                         "H4",
                                         ifelse(grepl("H8_", po4_sam),
                                                "H8",
                                                visit_id))))) %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H0_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H1_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H4_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H8_",
            replacement = "")


# Isolate rows with combined visit id
new_rows <- 
  all_phosphate2 %>% 
  filter(visit_id == "H0H1H4H8")

# Remove rows from dataframe and add mutated ones
all_phosphate2 <-
  all_phosphate2 %>% 
  anti_join(new_rows) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H0")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H1")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H4")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H8"))

# Get rid of blank wells
all_phosphate2 <- 
  all_phosphate2 %>% 
  rename(bromeliad_id = po4_sam) %>% 
  unite(column, po4_stds, bromeliad_id, sep = "-", remove = T) %>% 
  filter(column != "NA-NA") %>% 
  separate(column, into = c("po4_stds", "bromeliad_id"), sep = "-", remove = TRUE)

#replace "NA" with actual NA
all_phosphate2[all_phosphate2 == "NA"] = NA

# Check on structure of variables & make adjustments as necessary
str(all_phosphate2)
all_phosphate2 <- 
  all_phosphate2 %>% 
  mutate_at(vars(po4_stds), ~as.numeric(.)) %>% 
  mutate_at(vars(date, bromeliad_id, visit_id), ~as.character(.))

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
po4_split <- 
  all_phosphate2 %>%
  group_by(visit_id, plate) %>% 
  group_split()

# Create data frame with zeroed cal curve values for each plate  
plate_zeroes_po4 <- 
  po4_split %>% 
  map_df(plate_specific_zeroes_po4)

# Small dataframe with only plate zero values
summ_plate_zeroes_po4 <- 
  plate_zeroes_po4 %>% 
  group_by(visit_id, plate) %>% 
  summarise(mean_zero = mean(plate_zero))

# New column to sort by and nest dataset so it is compatible with the map function
nest_by_week_plate_po4 <- 
  plate_zeroes_po4 %>% 
  group_by(visit_id, plate) %>% 
  nest() 

# Apply function to each plate to get data frame with all of your slopes and intercepts
phosphate_mod_parms <- 
  nest_by_week_plate_po4 %>% 
  mutate(fit = map(data, po4_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(visit_id, plate, tidied) %>% 
  unnest(cols = c(tidied)) %>% 
  dplyr::select(visit_id, plate, term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -visit_id, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = po4_stds)

# Compare experimental values to expected ones
# Add plate zeroes to parameter table
phosphate_mod_parms1  <- 
  left_join(phosphate_mod_parms, summ_plate_zeroes_po4)

# Getting ready to plop expected parameter values into this table
pho_parms <- 
  nut_parms %>% 
  filter(nuts == "po4")

# Make vector for expected intercept to merge with modeled parameters table
pho_exp_int <- 
  rep(pho_parms$b, length = length(phosphate_mod_parms$slope))

# Make vector for expected slope to merge with modeled parameters table
pho_exp_slope <- 
  rep(pho_parms$m, length = length(phosphate_mod_parms$slope))

# Add expected values to the parm table
phosphate_mod_parms2 <- 
  cbind(phosphate_mod_parms1, pho_exp_int) %>% 
  rename(pho_exp_int = ...6)
phosphate_mod_parms3 <- 
  cbind(phosphate_mod_parms2, pho_exp_slope) %>% 
  rename(pho_exp_slope = ...7)

# Add columns for degree of divergence from expected vlaues
phosphate_mod_parms4 <- 
  phosphate_mod_parms3 %>% 
  mutate(int_dist = abs(pho_exp_int - intercept),
         slope_dist = abs(pho_exp_slope - slope))

# Replace intercept value with plate zero value if intercept is too far from expected
phosphate_mod_parms5 <- 
  phosphate_mod_parms4 %>% 
  mutate(intercept = ifelse(int_dist < 0.4*pho_exp_int, intercept, mean_zero))

# Recalculate int_dist with new intercept value
phosphate_mod_parms6 <- 
  phosphate_mod_parms5 %>% 
  mutate(int_dist = pho_exp_int - intercept)

# If zero values are still too far from expected, throw out run
pho_dist_int <- 
  phosphate_mod_parms6 %>% 
  filter(int_dist > 0.4*pho_exp_int)

# Look at slopes
pho_dist_slope <- 
  phosphate_mod_parms6 %>% 
  filter(slope_dist > 0.4*pho_exp_slope)

# Predict concentrations with new parm values
pho_sams <- 
  all_phosphate2 %>% 
  group_by(visit_id, plate) %>% 
  filter(!is.na(bromeliad_id)) %>% 
  dplyr::select(-po4_stds)

# Calculate concentrations!
pho_sams1 <- 
  left_join(pho_sams, phosphate_mod_parms6) %>% 
  dplyr::select(-pho_exp_int, -pho_exp_slope, -int_dist, -slope_dist) %>% 
  mutate(po4_con = (absorbance - intercept) / slope) 

# Calculate LOD
#LOD = 3*sd(slope)/slope
pho_sd_m <- 
  phosphate_mod_parms6 %>% 
  group_by(visit_id) %>% 
  summarise(sd_m = sd(slope))

pho_mean_m <- 
  phosphate_mod_parms6 %>% 
  group_by(visit_id) %>%
  summarise(mean_m = mean(slope)) %>% 
  left_join(pho_sd_m) %>% 
  mutate(po4_lod = 3*(sd_m)/(mean_m))

# Make negative concentrations and and positives below the LOD zeroes
po4_values <- 
  pho_sams1 %>% 
  left_join(pho_mean_m %>% 
              dplyr::select(visit_id, po4_lod)) %>% 
  mutate(po4_con = ifelse(po4_con < 0 | po4_con < po4_lod, 0, po4_con)) %>% 
  pivot_wider(values_from = po4_con, names_from = plate) %>% 
  rename(po4_a = A,
         po4_b = B)


# Ammonium ----------------------------------------------------------------
# Add in replicate information
all_ammonium2 <- 
  all_ammonium %>% 
  gather(key = round, value = absorbance, nh4_a, nh4_b) %>% 
  mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  left_join(visits %>% 
              dplyr::select(visit_id, date))%>% 
  mutate(visit_id = ifelse(grepl("H0_", nh4_sam),
                           "H0", 
                           ifelse(grepl("H1_", nh4_sam),
                                  "H1",
                                  ifelse(grepl("H4_", nh4_sam),
                                         "H4",
                                         ifelse(grepl("H8_", nh4_sam),
                                                "H8",
                                                visit_id))))) %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H0_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H1_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H4_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H8_",
            replacement = "")

# Isolate rows with combined visit id
new_rows <- 
  all_ammonium2 %>% 
  filter(visit_id == "H0H1H4H8")

# Remove rows from dataframe and add mutated ones
all_ammonium2 <-
  all_ammonium2 %>% 
  anti_join(new_rows) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H0")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H1")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H4")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H8"))

# Get rid of blank wells
all_ammonium2  <- 
  all_ammonium2 %>% 
  rename(bromeliad_id = nh4_sam) %>% 
  unite(column, nh4_stds, bromeliad_id, sep = "-", remove = T) %>% 
  filter(column != "NA-NA") %>% 
  separate(column, into = c("nh4_stds", "bromeliad_id"), sep = "-", remove = TRUE)

# Replace "NA" with actual NA
all_ammonium2[all_ammonium2 =="NA"] = NA

# C heck on structure of variables & make adjustments as necessary
str(all_ammonium2)
all_ammonium2 <- 
  all_ammonium2 %>% 
  mutate_at(vars(nh4_stds), ~as.numeric(.)) %>% 
  mutate_at(vars(date, bromeliad_id, visit_id), ~as.character(.)) 

# Some quick plotting to see how my calibration curves look
all_ammonium2 %>% 
  filter(!is.na(nh4_stds)) %>% 
  ggplot(aes(x = nh4_stds, y = absorbance)) +
  geom_point() +
  facet_grid(visit_id ~ plate)

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
nh4_split <- 
  all_ammonium2 %>%
  group_by(visit_id, plate) %>% 
  group_split()

# Make data frame with zeroed cal curve values for each plate  
plate_zeroes_nh4 <- 
  nh4_split %>% 
  map_df(plate_specific_zeroes_nh4)

# Get plate values
summ_plate_zeroes_nh4 <- 
  plate_zeroes_nh4 %>% 
  group_by(visit_id, plate) %>% 
  summarise(mean_zero = mean(plate_zero))


# New column to sort by and nest dataset so it is compatible with the map function
nest_by_date_plate_nh4 <- 
  plate_zeroes_nh4 %>% 
  group_by(visit_id, plate) %>% 
  nest() 


# Apply function to each plate to get data frame with each slope and intercept
ammonium_mod_parms <- 
  nest_by_date_plate_nh4 %>% 
  mutate(fit = map(data, nh4_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(visit_id, plate, tidied) %>% 
  unnest(cols = c(tidied)) %>% 
  dplyr::select(visit_id, plate, term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -visit_id, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = nh4_stds)

# Compare experimental parm values to expected ones 
# Add plate zeroes to parameter table
ammonium_mod_parms  <- 
  left_join(ammonium_mod_parms, summ_plate_zeroes_nh4)

# Getting ready to plop expected parameter values into this table
amm_parms <- 
  nut_parms %>% 
  filter(nuts == "nh4")

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
  filter(int_dist > 0.4*amm_exp_int) #looks good -- no data in this set

# Look at slopes
nh4_dist_slope <- 
  ammonium_mod_parms %>% 
  filter(slope_dist > 0.4*amm_exp_slope) 

# Predict concentrations with new parm values 
amm_split <- 
  all_ammonium2 %>% 
  group_by(visit_id, plate) %>% 
  filter(!is.na(bromeliad_id)) %>% 
  dplyr::select(-nh4_stds)

# Calculate concentrations!
amm_split <- 
  left_join(amm_split, ammonium_mod_parms) %>% 
  dplyr::select(-amm_exp_int, -amm_exp_slope, -int_dist, -slope_dist, -Wells) %>% 
  mutate(nh4_con = (absorbance - intercept) / slope)

# Calculate LODs
# LOD = 3*sd(slope)/slope
amm_sd_m <- 
  ammonium_mod_parms %>% 
  group_by(visit_id) %>% 
  summarise(sd_m = sd(slope))

amm_mean_m <- 
  ammonium_mod_parms %>% 
  group_by(visit_id) %>%
  summarise(mean_m = mean(slope)) %>% 
  left_join(amm_sd_m) %>% 
  mutate(nh4_lod = 3*(sd_m)/(mean_m))


# Negative concentrations = 0, concentrations under the detection limit = 0
nh4_values <- 
  amm_split %>% 
  left_join(amm_mean_m %>% 
              dplyr::select(visit_id, nh4_lod)) %>% 
  mutate(nh4_con = ifelse(nh4_con < 0 | nh4_con < nh4_lod, 
                          0, nh4_con)) %>% 
  pivot_wider(values_from = nh4_con, names_from = plate) %>% 
  rename(nh4_a = A,
         nh4_b = B)



# Bacteria ----------------------------------------------------------------
# Add in replicate information and actual concentration of spheres in standard
all_bacteria2 <- 
  all_bacteria %>% 
  gather(key = round, value = absorbance, bact_a, bact_b) %>% 
  mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  dplyr::select(-round) %>% 
  left_join(visits %>% 
              dplyr::select(visit_id, date, spheres_per_L)) %>% 
  mutate(bact_stds = ifelse(!is.na(bact_stds), 
                            bact_stds * spheres_per_L, bact_stds)) %>% 
  dplyr::select(-spheres_per_L)%>% 
  mutate(visit_id = ifelse(grepl("H0_", bact_sam),
                           "H0", 
                           ifelse(grepl("H1_", bact_sam),
                                  "H1",
                                  ifelse(grepl("H4_", bact_sam),
                                         "H4",
                                         ifelse(grepl("H8_", bact_sam),
                                                "H8",
                                                visit_id))))) %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H0_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H1_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H4_",
            replacement = "") %>% 
  mutate_if(is.character,
            stringr::str_replace_all, 
            pattern = "H8_",
            replacement = "")

# Isolate rows with combined visit id
new_rows <- 
  all_bacteria2 %>% 
  filter(visit_id == "H0H1H4H8")

# Remove rows from dataframe and add mutated ones
all_bacteria2 <-
  all_bacteria2 %>% 
  anti_join(new_rows) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H0")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H1")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H4")) %>% 
  bind_rows(new_rows %>% 
              mutate(visit_id = "H8"))

# Get rid of blank wells
all_bacteria2  <- 
  all_bacteria2 %>% 
  rename(bromeliad_id = bact_sam) %>% 
  unite(column, bact_stds, bromeliad_id, sep = "-", remove = T) %>% 
  filter(column != "NA-NA") %>% 
  separate(column, into = c("bact_stds", "bromeliad_id"), sep = "-", remove = TRUE)

# Replace "NA" with actual NA
all_bacteria2[all_bacteria2 =="NA"] = NA

# Check on structure of variables & make adjustments as necessary
str(all_bacteria2)
all_bacteria2 <- 
  all_bacteria2 %>% 
  mutate_at(vars(bact_stds), ~as.numeric(.)) %>% 
  mutate_at(vars(date, bromeliad_id, visit_id), ~as.character(.)) 

# Some quick plotting to see how my calibration curves look
all_bacteria2 %>% 
  filter(!is.na(bact_stds)) %>% 
  ggplot(aes(x = bact_stds, y = absorbance)) +
  geom_point() +
  facet_grid(visit_id ~ plate)

# Split by date and plate, so that each plate's samples are fit to the proper cal curve
bact_split <- 
  all_bacteria2 %>%
  group_by(visit_id, plate) %>% 
  group_split()

# Make data frame with zeroed cal curve values for each plate  
plate_zeroes_bact <- 
  bact_split %>% 
  map_df(plate_specific_zeroes_bact)

# Get plate values
summ_plate_zeroes_bact <- 
  plate_zeroes_bact %>% 
  group_by(visit_id, plate) %>% 
  summarise(mean_zero = mean(plate_zero))


# New column to sort by and nest dataset so it is compatible with the map function
nest_by_date_plate_bact <- 
  plate_zeroes_bact %>% 
  group_by(visit_id, plate) %>% 
  nest() 


# Apply function to each plate to get data frame with each slope and intercept
bacteria_mod_parms <- 
  nest_by_date_plate_bact %>% 
  mutate(fit = map(data, bact_cal_curve),
         tidied = map(fit, tidy)) %>% 
  dplyr::select(visit_id, plate, tidied) %>% 
  unnest(cols = c(tidied)) %>% 
  dplyr::select(visit_id, plate, term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -visit_id, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = bact_stds)


# Predict concentrations with new parm values 
bact_split <- 
  all_bacteria2 %>% 
  group_by(visit_id, plate) %>% 
  filter(!is.na(bromeliad_id)) %>% 
  dplyr::select(-bact_stds)

# Calculate concentrations!
bact_split <- 
  left_join(bact_split, bacteria_mod_parms) %>% 
    mutate(bact_con = (absorbance - intercept) / slope)



# Negative concentrations = 0, concentrations under the detection limit = 0
## First visit NA because different, uncomparable method
bact_values <- 
  bact_split %>% 
  mutate(bact_con = ifelse(visit_id ==1, NA,
                           ifelse(bact_con < 0 & visit_id > 1,
                                  0, bact_con))) %>%
  pivot_wider(values_from = bact_con, names_from = plate) %>% 
  rename(bact_a = A,
         bact_b = B)


# Convert chlorophyll values ----------------------------------------------
chlorophyll <- 
  chlorophyll %>% 
  ## Remove blank values
  left_join(chlorophyll %>% 
              filter(bromeliad_id == "acetone") %>% 
              dplyr::select(visit_id, chlorophyll_raw) %>% 
              group_by(visit_id) %>% 
              summarise_all(mean) %>% 
              rename(blanks = chlorophyll_raw)) %>% 
  mutate(chlorophyll_raw = chlorophyll_raw - blanks) %>% 
  dplyr::select(-blanks) %>% 
  mutate(chlorophyll_ugL = chlorophyll_raw * (solvent_added / v_filtered))

  
  
# Combine and save data ---------------------------------------------------
weekly_measurements <- 
  weekly_measurements_init %>%
  mutate(visit_id = as.character(visit_id),
         bromeliad_id = as.character(bromeliad_id)) %>% 
  left_join(bromeliads) %>% 
  left_join(subsidies %>% 
              mutate(bromeliad_id = as.character(bromeliad_id)) %>% 
              dplyr::select(bromeliad_id, litter_1_mg, feces_1_mg, litter_2_mg, feces_2_mg, litter_3_mg, feces_3_mg)) %>% 
  left_join(chlorophyll %>% 
              dplyr::select(visit_id, bromeliad_id, chlorophyll_ugL) %>% 
              mutate(visit_id = as.character(visit_id),
                     bromeliad_id = as.character(bromeliad_id))) %>% 
  left_join(po4_values %>% 
              dplyr::select(visit_id, bromeliad_id, po4_a) %>% 
              filter(!is.na(po4_a)) %>% 
              unique()) %>% 
  left_join(po4_values %>% 
              dplyr::select(visit_id, bromeliad_id, po4_b) %>% 
              filter(!is.na(po4_b)) %>%
              unique()) %>%
  left_join(nh4_values %>% 
              dplyr::select(visit_id, bromeliad_id, nh4_a) %>% 
              filter(!is.na(nh4_a)) %>% 
              unique()) %>% 
  left_join(nh4_values %>% 
              dplyr::select(visit_id, bromeliad_id, nh4_b) %>% 
              filter(!is.na(nh4_b)) %>%
              unique()) %>%
  left_join(no2_values %>% 
              dplyr::select(visit_id, bromeliad_id, no2_a) %>% 
              filter(!is.na(no2_a)) %>%
              unique()) %>% 
  left_join(no2_values %>% 
            dplyr::select(visit_id, bromeliad_id, no2_b) %>% 
            filter(!is.na(no2_b)) %>%
            unique()) %>% 
  left_join(no3_values %>% 
              dplyr::select(visit_id, bromeliad_id, no3_a) %>% 
              mutate(visit_id = as.character(visit_id),
                     bromeliad_id = as.character(bromeliad_id)) %>% 
              filter(!is.na(no3_a)) %>%
              unique()) %>% 
  left_join(no3_values %>% 
              dplyr::select(visit_id, bromeliad_id, no3_b) %>% 
              mutate(visit_id = as.character(visit_id),
                     bromeliad_id = as.character(bromeliad_id)) %>% 
              filter(!is.na(no3_b)) %>%
              unique()) %>% 
  left_join(bact_values %>% 
              dplyr::select(visit_id, bromeliad_id, bact_a) %>% 
              filter(!is.na(bact_a)) %>% 
              unique()) %>% 
  left_join(bact_values %>% 
              dplyr::select(visit_id, bromeliad_id, bact_b) %>% 
              filter(!is.na(bact_b)) %>%
              unique()) %>% 
  ## Average two values
  dplyr::rowwise() %>% 
  dplyr::mutate(po4 = mean(c(po4_a, po4_b), na.rm = T),
                no2 = mean(c(no2_a, no2_b), na.rm = T),
                nh4 = mean(c(nh4_a, nh4_b), na.rm = T),
                no3 = mean(c(no3_a, no3_b), na.rm = T),
                bact = mean(c(bact_a, bact_b), na.rm = T)) %>% 
  dplyr::select(-po4_a, -po4_b, -nh4_a, -nh4_b,
                -no2_a, -no2_b, -no3_a, -no3_b,
                -bact_a, -bact_b, -v_extracted) %>% 
  ## Add subsidy
  mutate(subsidy_1 = ifelse(feces_1_mg == 0, 
                            "litter_only", "litter_feces"),
         subsidy_2 = ifelse(feces_2_mg == 0, 
                            "litter_only", "litter_feces"),
         subsidy_3 = ifelse(feces_3_mg == 0, 
                            "litter_only", "litter_feces")) %>% 
  unite(subsidy,subsidy_1, subsidy_2, sep = "_", remove = F) %>% 
  dplyr::mutate(subsidy = ifelse(subsidy == "NA_NA", subsidy_3,
                                 subsidy)) %>% 
  dplyr::select(-litter_1_mg, -litter_2_mg, -litter_3_mg, -subsidy_2,
                -feces_1_mg, -feces_2_mg, -feces_3_mg, -subsidy_3) %>% 
  ## Add date
  dplyr::left_join(visits %>% 
                     dplyr::select(visit_id, date)) %>% 
  dplyr::mutate(date = lubridate::yday(lubridate::as_date(lubridate::ymd(date)))) 

write.csv(weekly_measurements,
          "microchannels/appdata/bromeliad_tax_exp.csv",
          row.names = F)

# Assertions ---------------------------------------------------------------------
# Check data format
str(weekly_measurements)

# Make sure all data are in OK range
weekly_measurements %>% 
  mutate(n_leaves = as.numeric(n_leaves),
         larvae = as.numeric(larvae),
         subsidy_id = as.numeric(subsidy_id),
         capacity_mL = as.numeric(capacity_mL)) %>% 
  dplyr::select(-vessel, -shading, -bromeliad_id, -visit_id) %>% 
  assert(within_bounds(0, Inf),
         everything()) %>% 
  verify(nrow(.) == 312)

# Make sure no typos in several columns
## 12 visits
length(unique(weekly_measurements$visit_id)) == 12
## 52 bromelias
length(unique(weekly_measurements$bromeliad_id)) == 52
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


# SessionInfo -------------------------------------------------------------
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 22621)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.1   car_3.1-0       carData_3.0-5   ggeffects_1.1.2 here_1.0.1     
# [6] forcats_0.5.2   stringr_1.4.1   dplyr_1.1.0     purrr_0.3.4     readr_2.1.3    
# [11] tidyr_1.2.0     tibble_3.1.8    ggplot2_3.4.0   tidyverse_1.3.2 lmerTest_3.1-3 
# [16] lme4_1.1-31     Matrix_1.4-1   
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.4.4          jsonlite_1.8.3      splines_4.2.1       modelr_0.1.9       
# [5] assertthat_0.2.1    stats4_4.2.1        tensorA_0.36.2      googlesheets4_1.0.1
# [9] cellranger_1.1.0    pbivnorm_0.6.0      numDeriv_2016.8-1.1 pillar_1.8.1       
# [13] backports_1.4.1     lattice_0.20-45     glue_1.6.2          rvest_1.0.3        
# [17] minqa_1.2.5         colorspace_2.0-3    pkgconfig_2.0.3     broom_1.0.0        
# [21] haven_2.5.0         corpcor_1.6.10      scales_1.2.1        tzdb_0.3.0         
# [25] cubature_2.0.4.5    googledrive_2.0.0   mgcv_1.8-42         generics_0.1.3     
# [29] ellipsis_0.3.2      withr_2.5.0         cli_3.5.0           mnormt_2.1.0       
# [33] crayon_1.5.2        magrittr_2.0.3      readxl_1.4.1        fs_1.5.2           
# [37] fansi_1.0.3         nlme_3.1-157        MASS_7.3-57         xml2_1.3.3         
# [41] tools_4.2.1         hms_1.1.2           gargle_1.2.0        lifecycle_1.0.3    
# [45] munsell_0.5.0       reprex_2.0.2        DHARMa_0.4.5        compiler_4.2.1     
# [49] rlang_1.0.6         grid_4.2.1          nloptr_2.0.3        rstudioapi_0.13    
# [53] lavaan_0.6-12       boot_1.3-28         gtable_0.3.1        abind_1.4-5        
# [57] DBI_1.1.3           R6_2.5.1            lubridate_1.8.0     utf8_1.2.2         
# [61] rprojroot_2.0.3     MCMCglmm_2.33       ape_5.6-2           stringi_1.7.8      
# [65] parallel_4.2.1      Rcpp_1.0.9          vctrs_0.5.2         dbplyr_2.2.1       
# [69] tidyselect_1.2.0    coda_0.19-4      

