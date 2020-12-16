#processing nutrients for the Nact experiment from TEAMM 2018. It is a multi-plate processing script. It requires that all of the datafiles are in plater-compatible format (see: https://cran.r-project.org/web/packages/plater/vignettes/plater-basics.html) and are stored in a single folder
#all of these concentrations are in ug/L


#load necessary packages
library(plater)
library(tidyverse)
library(stringr)
library(cowplot)
library(modelr)
library(lubridate)
library(broom)

#load in standard parameter values dataset
nut_parms <- read_csv(file = "data/Nact/nut_parms.csv")

#get the data ----------------------------------------------------
#read in the data -- here we are making a list of datafiles
nutrient_files <-  c(list.files("data/Nact/nuts", full.names = TRUE))

#sort files by nutrient
nitrate_files <- nutrient_files[grepl(".no", nutrient_files)]
phosphate_files <- nutrient_files[grepl(".po", nutrient_files)]
ammonium_files <- nutrient_files[grepl(".nh", nutrient_files)]

#create function to read in plate-shaped csv file, make into dataframe, rename columns
read_plate_function <- function(file_name) {
  data <- read_plate(
    file = file_name,            
    well_ids_column = "Wells"  
  )
}

#apply above function to each data file within nitrate_files --> output is one large data frame for each nutrient with data from all dates
all_nitrate <- nitrate_files %>% 
  map_df(read_plate_function)
all_phosphate <- phosphate_files %>% 
  map_df(read_plate_function)
all_ammonium <- ammonium_files %>% 
  map_df(read_plate_function)

#nitrate DONE 13 MARCH------------------------------------------------------------------

#### data tidying ####
str(all_nitrate)

#separate by plate -- each plate it a week's worth of data
nit_wk1 <- all_nitrate[1:96,]
nit_wk2 <- all_nitrate[97:192,]

#create week column for each new dataframe
nit_wk1 <- nit_wk1 %>% 
  mutate(week = "nit1")
nit_wk2 <- nit_wk2 %>% 
  mutate(week = "nit2")

#glue them back together
all_nitrate1 <- rbind(nit_wk1, nit_wk2)

#make data long form, add plate column, change round to just be no2 or no2no3
all_nitrate2 <- all_nitrate1 %>% 
  gather(key = round, value = absorbance, no2_a, no2_b, no2no3_a, no2no3_b) %>% 
  mutate(replicate = if_else(grepl("a$", round), "A", "B")) %>% 
  mutate(Round = if_else(grepl("no3", round), "no2no3", "no2")) %>%
  select(-round) %>% 
  rename(round = Round) %>% 
  spread(key = round, value = absorbance)

#add column for net absorbance due to nitrate
all_nitrate3 <- all_nitrate2 %>% 
  mutate(no3 = no2no3 - no2)

#check on structure of variables & make adjustments as necessary
str(all_nitrate3)
all_nitrate3 <- all_nitrate3 %>% 
  mutate_at(vars(replicate, date, no2no3_sam, week), ~as.factor(.)) %>% 
  filter(no3 < 0.5) %>% 
  select(-Wells)

#some quick plotting to see how my calcurves look -- zeroes look wonky for nit2
all_nitrate3 %>% 
  filter(!is.na(no3_stds)) %>% 
  ggplot(aes(x = no3_stds, y = no3)) +
  geom_point() +
  facet_grid(week ~ replicate)

#split by date and plate, so that each plate's samples are fit to the proper cal curve
nit_split <- all_nitrate3 %>%
  group_by(week, replicate) %>% 
  group_split()

#### get cal curve parameters for each plate ####
#zero each plate according to its own blank values
plate_specific_zeroes_no3 <- function(df){
  standards <- df %>% 
    filter(!is.na(no3_stds))  #isolate no3_standards
  blank_value_no3 <- standards %>% #get blank value
    filter(no3_stds < 0.1) %>% 
    summarise(mean(no3)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_no3)
}

#create data frame with zeroed cal curve values for each plate  
full_plate_zeroes_no3 <- nit_split %>% 
  map_df(plate_specific_zeroes_no3)
  
#small df with only plate zero values
summ_plate_zeroes_no3 <- full_plate_zeroes_no3 %>% 
  group_by(week, replicate) %>% 
  summarise(mean_zero = mean(plate_zero))

#new column to sort by and nest dataset so it is compatible with the map function
nest_by_week_rep_no3 <- plate_zeroes_no3 %>% 
  group_by(week, replicate) %>% 
  nest() 

#write function to make cal curves
no3_cal_curve <- function(df) {
  lm(no3 ~ no3_stds, data = df)
}

#apply cal curve function to each plate to get data frame with all slopes and intercepts
nitrate_mod_parms <- nest_by_week_rep_no3 %>% 
  mutate(fit = map(data, no3_cal_curve),
         tidied = map(fit, tidy)) %>% 
  select(week, replicate, tidied) %>% 
  unnest() %>% 
  select(week, replicate,  term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -week, -replicate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = no3_stds)


#### compare exp parm values to exp ones ####
nitrate_mod_parms  <- left_join(nitrate_mod_parms, summ_plate_zeroes_no3)

#getting ready to plop expected parameter values into this table
nit_parms <- nut_parms %>% 
  filter(nuts == "no3")

#make vector for expected intercept to merge with modeled parameters table
nit_exp_int <- rep(nit_parms$b, length = length(nitrate_mod_parms$intercept))

#make vector for expected slope to merge with modeled parameters table
nit_exp_slope <- rep(nit_parms$m, length = length(nitrate_mod_parms$intercept))

#adding expected values to the parm table
nitrate_mod_parms <- cbind(nitrate_mod_parms, nit_exp_int)
nitrate_mod_parms <- cbind(nitrate_mod_parms, nit_exp_slope)

#add columns for degree of divergence from expected vlaues
nitrate_mod_parms <- nitrate_mod_parms %>% 
  mutate(int_dist = abs(nit_exp_int - intercept),
         slope_dist = abs(nit_exp_slope - slope))

#replace intercept value with plate zero value if intercept is too far from expected
nitrate_mod_parms <- nitrate_mod_parms %>% 
  mutate(intercept = ifelse(int_dist < 0.3*nit_exp_int, intercept, mean_zero))

#recalculate int_dist with new intercept values
nitrate_mod_parms <- nitrate_mod_parms %>% 
  mutate(int_dist = nit_exp_int - intercept)

#if zero values are still too far from expected, throw out run -- good to go
nit_dist_int <- nitrate_mod_parms %>% 
  filter(int_dist > 0.3*nit_exp_int)

#look at slopes -- slopes look good to go
nit_dist_slope <- nitrate_mod_parms %>% 
  filter(abs(slope_dist) > 0.4*nit_exp_slope)

#LOD = 3*sd(slope)/slope
#nitrate
nit_sd_m <- 
  nitrate_mod_parms %>% 
  group_by(visit_id) %>% 
  summarise(sd_m = sd(slope)) 

nit_mean_m <- 
  nitrate_mod_parms %>% 
  group_by(visit_id) %>%
  summarise(mean_m = mean(slope)) 

nit_LOD <- 3*(nit_sd_m)/(nit_mean_m)


nitrate_mod_parms <- 
  nitrate_mod_parms %>% 
  mutate(no3_lod = nit_LOD)

#### predict concentrations with new parm values ####
nit_sams <- all_nitrate3 %>% 
  group_by(week, replicate) %>% 
  filter(!is.na(no2no3_sam)) %>% 
  select(-no3_stds, -no2_stds)

#use lm parameters for each plate to estimate concentration
nit_sams1 <- left_join(nit_sams, nitrate_mod_parms) %>% 
  select(-nit_exp_int, -nit_exp_slope, -int_dist, -slope_dist) %>% 
  mutate(`[no3]` = (no3 - intercept) / slope) %>% 
  mutate(`[no3]` = ifelse(`[no3]` < 0, 0, `[no3]`)) %>% 
  mutate(`[no3]` = ifelse(`[no3]` < nit_LOD, 0, `[no3]`))#figure out what to do with negatives

#new data frame with tank-date averages
nit_sams_avg <- nit_sams1 %>%
  group_by(no2no3_sam, date) %>% 
  summarise(mean_no3 = mean(`[no3]`))

#export nitrate concentrations to csv
#write_csv(nit_sams_avg, "data/Nact/nuts/nitrate_levels_processed.csv")

#quick plotting -- this is not working properly
nit_sams_avg %>% 
  group_by(no2no3_sam) %>% 
  ggplot(aes(x = as.numeric(as.character(date)), y = mean_no3)) + 
  geom_point() + 
  geom_line()


#phosphate DONE 14 MARCH -------------------------------------------------------------

#### data tidying ####

pho_1 <- all_phosphate[1:96,]
pho_2 <- all_phosphate[97:192,]

#create unique ID column for each new dataframe
pho_1 <- pho_1 %>% 
  mutate(week = "1")
pho_2 <- pho_2 %>% 
  mutate(week = "2")

#glue them back together
all_phosphate1 <- rbind(pho_1, pho_2)
all_phosphate1 <- all_phosphate1 %>% 
  gather(key = round, value = absorbance, po4_a, po4_b) %>% 
  mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  select(-round)

#get rid of blank wells
all_phosphate2 <- all_phosphate1 %>% 
  unite(column, po4_stds, po4_sam, sep = "-", remove = T) %>% 
  filter(column != "NA-NA") %>% 
  separate(column, into = c("po4_stds", "po4_sam"), sep = "-", remove = TRUE)

#replace "NA" with actual NA
all_phosphate2[all_phosphate2 =="NA"] = NA

#check on structure of variables & make adjustments as necessary
str(all_phosphate2)
all_phosphate2 <- all_phosphate2 %>% 
  mutate_at(vars(po4_stds), ~as.numeric(.)) %>% 
  mutate_at(vars(date, po4_sam, week), ~as.factor(.))

#some quick plotting to see how my calcurves look -- zeroes look wonky for nit2
all_phosphate2 %>% 
  filter(!is.na(po4_stds)) %>% 
  ggplot(aes(x = po4_stds, y = absorbance)) +
  geom_point() +
  facet_grid(week ~ plate) #week 1B slope may be low


#split by date and plate, so that each plate's samples are fit to the proper cal curve
po4_split <- all_phosphate2 %>%
  group_by(week, plate) %>% 
  group_split()

#### get cal curve parameters for each plate ####
#zero each plate according to its own blank values
plate_specific_zeroes_po4 <- function(df){
  standards <- df %>% 
    filter(!is.na(po4_stds))  #isolate no3_standards
  blank_value_po4 <- standards %>% #get blank value
    filter(po4_stds == 0.0) %>% 
    summarise(mean(absorbance)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_po4)
}

#create data frame with zeroed cal curve values for each plate  
plate_zeroes_po4 <- po4_split %>% 
  map_df(plate_specific_zeroes_po4)

#small df with only plate zero values
summ_plate_zeroes_po4 <- plate_zeroes_po4 %>% 
  group_by(week, plate) %>% 
  summarise(mean_zero = mean(plate_zero))

#new column to sort by and nest dataset so it is compatible with the map function
nest_by_week_plate_po4 <- plate_zeroes_po4 %>% 
  group_by(week, plate) %>% 
  nest() 

#write function to make cal curves
po4_cal_curve <- function(df) {
  lm(absorbance ~ po4_stds, data = df)
}

#apply function to each plate to get data frame with all of your slopes and intercepts
phosphate_mod_parms <- nest_by_week_plate_po4 %>% 
  mutate(fit = map(data, po4_cal_curve),
         tidied = map(fit, tidy)) %>% 
  select(week, plate, tidied) %>% 
  unnest() %>% 
  select(week, plate, term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  gather(key = term, value = estimate, -week, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = po4_stds)

#### compare exp parm values to exp ones ####
#add plate zeroes to parameter table
phosphate_mod_parms1  <- left_join(phosphate_mod_parms, summ_plate_zeroes_po4)

#getting ready to plop expected parameter values into this table
pho_parms <- nut_parms %>% 
  filter(nuts == "po4")

#make vector for expected intercept to merge with modeled parameters table
pho_exp_int <- rep(pho_parms$b, length = length(phosphate_mod_parms$slope))

#make vector for expected slope to merge with modeled parameters table
pho_exp_slope <- rep(pho_parms$m, length = length(phosphate_mod_parms$slope))

#adding expected values to the parm table
phosphate_mod_parms2 <- cbind(phosphate_mod_parms1, pho_exp_int)
phosphate_mod_parms3 <- cbind(phosphate_mod_parms2, pho_exp_slope)

#add columns for degree of divergence from expected vlaues
phosphate_mod_parms4 <- phosphate_mod_parms3 %>% 
  mutate(int_dist = abs(pho_exp_int - intercept),
         slope_dist = abs(pho_exp_slope - slope))

#replace intercept value with plate zero value if intercept is too far from expected
phosphate_mod_parms5 <- phosphate_mod_parms4 %>% 
  mutate(intercept = ifelse(int_dist < 0.4*pho_exp_int, intercept, mean_zero))

#recalculate int_dist with new intercept value
phosphate_mod_parms6 <- phosphate_mod_parms5 %>% 
  mutate(int_dist = pho_exp_int - intercept)

#if zero values are still too far from expected, throw out run
pho_dist_int <- phosphate_mod_parms6 %>% 
  filter(int_dist > 0.4*pho_exp_int) #looks good -- no data in this set

#look at slopes
pho_dist_slope <- phosphate_mod_parms6 %>% 
  filter(slope_dist > 0.4*pho_exp_slope) #go back and look at these dates in your lab notebook -- these are all suspiciously close to 2x slope, so could have used a standard chemical with two phosphate molecules by accident #STILL NEED TO VERIFY

#if we assume this is the mistake, let's correct that and try again
phosphate_mod_parms7 <- phosphate_mod_parms6 %>% 
  mutate(slope = ifelse(slope < 0.015, slope, .5*slope))


#predict concentrations with new parm values --------------------------
pho_sams <- all_phosphate2 %>% 
  group_by(week, plate) %>% 
  filter(!is.na(po4_sam)) %>% 
  select(-po4_stds)

#calculate concentrations!
pho_sams1 <- left_join(pho_sams, phosphate_mod_parms7) %>% 
  select(-pho_exp_int, -pho_exp_slope, -int_dist, -slope_dist) %>% 
  mutate(`[po4]` = (absorbance - intercept) / slope) #figure out what to do with negatives


#calculate LOD
#LOD = 3*sd(slope)/slope
#phosphate
pho_sd_m <- phosphate_mod_parms7 %>% 
  summarise(sd_m = sd(slope)) %>% 
  unlist()

pho_mean_m <- phosphate_mod_parms7 %>% 
  summarise(mean_m = mean(slope)) %>% 
  unlist()

pho_LOD <- 3*(pho_sd_m)/(pho_mean_m)

#Make negative concentrations and and positives below the LOD zeroes
pho_sams2 <- pho_sams1 %>% 
  mutate(`[po4]` = ifelse(`[po4]` < 0, 0, `[po4]`),
         `[po4]` = ifelse(`[po4]` < pho_LOD, 0, `[po4]`))

#average phosphate for each date,tank pair
pho_sam_avg <- pho_sams2 %>% 
  group_by(po4_sam, date) %>% 
  summarise(mean_pho = mean(`[po4]`))
#ugh, these look NOT GOOD

#write to csv
#write_csv(pho_sam_avg, "data/Nact/nuts/phosphate_levels_processed.csv")

#quick plot
pho_sam_avg %>% 
  group_by(po4_sam) %>% 
  ggplot(aes(x = as.numeric(as.character(date)), y = mean_pho, colour = as.numeric(po4_sam))) +
  geom_point() #week 2 , middle metacommunities have some issues.
#need to re-run week 2 nuts

#


#ammonium ----------------------------------------------------------
#data tidying
amm_1 <- all_ammonium[1:96,]
amm_2 <- all_ammonium[97:192,]

#create unique ID column for each new dataframe
amm_1 <- amm_1 %>% 
  mutate(week = "1")
amm_2 <- amm_2 %>% 
  mutate(week = "2")

#glue them back together
all_ammonium1 <- rbind(amm_1, amm_2)

#add in replicate information
all_ammonium2 <- all_ammonium1 %>% 
  gather(key = round, value = absorbance, nh4_a, nh4_b) %>% 
  mutate(plate = if_else(grepl("a$", round), "A", "B")) %>% 
  select(-round)

#get rid of blank wells
all_ammonium3 <- all_ammonium2 %>% 
  unite(column, nh4_stds, nh4_sam, sep = "-", remove = T) %>% 
  filter(column != "NA-NA") %>% 
  separate(column, into = c("nh4_stds", "nh4_sam"), sep = "-", remove = TRUE)

#replace "NA" with actual NA
all_ammonium3[all_ammonium3 =="NA"] = NA

#check on structure of variables & make adjustments as necessary
str(all_ammonium3)
all_ammonium3 <- all_ammonium3 %>% 
  mutate_at(vars(nh4_stds), ~as.numeric(.)) %>% 
  mutate_at(vars(date, nh4_sam, plate, week), ~as.factor(.)) 

#some quick plotting to see how my calcurves look
all_ammonium3 %>% 
  filter(!is.na(nh4_stds)) %>% 
  ggplot(aes(x = nh4_stds, y = absorbance)) +
  geom_point() +
  facet_grid(week ~ plate)

#definitely an outlier in week 1 plate A-- throw that out
all_ammonium4 <- all_ammonium3 %>% 
  filter(absorbance < 0.2)

#split by date and plate, so that each plate's samples are fit to the proper cal curve
nh4_split <- all_ammonium4 %>%
  group_by(week, plate) %>% 
  group_split()

#get cal curve parameters for each plate --------------------------
#zero each plate according to its own blank values
plate_specific_zeroes_nh4 <- function(df){
  standards <- df %>% 
    filter(!is.na(nh4_stds))  #isolate standards
  blank_value_nh4 <- standards %>% #get blank value
    filter(nh4_stds == 0.0) %>% 
    summarise(mean(absorbance)) %>% 
    unlist()
  standards <- standards %>% 
    mutate(plate_zero = blank_value_nh4)
}

#create data frame with zeroed cal curve values for each plate  
plate_zeroes_nh4 <- nh4_split %>% 
  map_df(plate_specific_zeroes_nh4)

#get plate values
summ_plate_zeroes_nh4 <- plate_zeroes_nh4 %>% 
  group_by(week, plate) %>% 
  summarise(mean_zero = mean(plate_zero))


#new column to sort by and nest dataset so it is compatible with the map function
nest_by_date_plate_nh4 <- plate_zeroes_nh4 %>% 
  group_by(week, plate) %>% 
  nest() 

#write function to make cal curves
nh4_cal_curve <- function(df) {
  lm(absorbance ~ nh4_stds, data = df)
}

#apply function to each plate to get data frame with each slope and intercept
ammonium_mod_parms <- nest_by_date_plate_nh4 %>% 
  mutate(fit = map(data, nh4_cal_curve),
         tidied = map(fit, tidy)) %>% 
  select(week, plate, tidied) %>% 
  unnest() %>% 
  select(week, plate, term, estimate) %>% 
  spread(key = term, value = estimate) %>% 
  select(week, plate, `(Intercept)`, nh4_stds) %>% 
  gather(key = term, value = estimate, -week, -plate) %>% 
  filter(!is.na(estimate)) %>% 
  spread(key= term, value = estimate) %>% 
  rename(intercept = `(Intercept)`,
         slope = nh4_stds)

# compare exp parm values to exp ones --------------------------------
#add plate zeroes to parameter table
ammonium_mod_parms  <- left_join(ammonium_mod_parms, summ_plate_zeroes_nh4)

#getting ready to plop expected parameter values into this table
amm_parms <- nut_parms %>% 
  filter(nuts == "nh4")

#make vector for expected intercept to merge with modeled parameters table
amm_exp_int <- rep(amm_parms$b, length = length(ammonium_mod_parms$slope))

#make vector for expected slope to merge with modeled parameters table
amm_exp_slope <- rep(amm_parms$m, length = length(ammonium_mod_parms$slope))

#adding expected values to the parm table
ammonium_mod_parms <- cbind(ammonium_mod_parms, amm_exp_int)
ammonium_mod_parms <- cbind(ammonium_mod_parms, amm_exp_slope)

#add columns for degree of divergence from expected vlaues
ammonium_mod_parms <- ammonium_mod_parms %>% 
  mutate(int_dist = abs(amm_exp_int - intercept),
         slope_dist = abs(amm_exp_slope - slope))

#replace intercept value with plate zero value if intercept is too far from expected
ammonium_mod_parms <- ammonium_mod_parms %>% 
  mutate(intercept = ifelse(int_dist < 0.4*amm_exp_int, intercept, mean_zero))

#recalculate int_dist with new intercept value
ammonium_mod_parms <- ammonium_mod_parms %>% 
  mutate(int_dist = amm_exp_int - intercept)

#if zero values are still too far from expected, throw out run
nh4_dist_int <- ammonium_mod_parms %>% 
  filter(int_dist > 0.4*amm_exp_int) #looks good -- no data in this set

#look at slopes
#bad slopes --> go back and look at R-squared values (is there an outlier point in cal curve to delete)
nh4_dist_slope <- ammonium_mod_parms %>% 
  filter(slope_dist > 0.4*amm_exp_slope) #will need to look at/re-do week 2 plate A -- slope is low

#predict concentrations with new parm values ------------------------------
amm_split <- all_ammonium4 %>% 
  group_by(week, plate) %>% 
  filter(!is.na(nh4_sam)) %>% 
  select(-nh4_stds)

#calculate concentrations!
amm_split <- left_join(amm_split, ammonium_mod_parms) %>% 
  select(-amm_exp_int, -amm_exp_slope, -int_dist, -slope_dist) %>% 
  mutate(`[nh4]` = (absorbance - intercept) / slope)

#calculate LODs
#LOD = 3*sd(slope)/slope
#ammonium
amm_sd_m <- ammonium_mod_parms %>% 
  summarise(sd_m = sd(slope)) %>% 
  unlist()
amm_mean_m <- ammonium_mod_parms %>% 
  summarise(mean_m = mean(slope)) %>% 
  unlist()

amm_LOD <- 3*(amm_sd_m)/(amm_mean_m)

#negative concentrations = 0, concentrations under the dectection limit = 0
amm_split1 <- amm_split %>% 
  mutate(`[nh4]` = ifelse(`[nh4]` < 0, 0, `[nh4]`),
         `[nh4]` = ifelse(`[nh4]` < amm_LOD, 0, `[nh4]`))

#take averages for each tank-date
amm_sam_avg <- amm_split1 %>% 
  group_by(nh4_sam, date) %>% 
  summarise(mean_amm = mean(`[nh4]`))

#write results to csv
#write_csv(amm_sam_avg, "data/Nact/nuts/ammonium_levels_processed.csv")

# make whole NACT file ---------------------------------------------------
amm <- read_csv(file = "data/Nact/nuts/ammonium_levels_processed.csv")
nit <- read_csv(file = "data/Nact/nuts/nitrate_levels_processed.csv")
pho <- read_csv(file = "data/Nact/nuts/phosphate_levels_processed.csv")

#change each sample column to "tank"
amm1 <- amm %>% 
  rename(tank = nh4_sam)
nit1 <- nit %>% 
  rename(tank = no2no3_sam)
pho1 <- pho %>% 
  rename(tank = po4_sam)

#combine together
nact_nuts <- left_join(nit1, pho1)
nact_nuts <- left_join(nact_nuts, amm1)

#make the date look like a date so it is easier to merge with other data
nact_nuts1 <- nact_nuts %>% 
  separate(date, into = c("year", "month"), sep = 4) %>% 
  separate(month, into = c("month", "day"), sep = 2) %>% 
  unite(date, year, month, day, sep = "-") %>% 
  mutate(date = ymd(date)) %>% 
  rename(nitrate = mean_no3,
         ammonium = mean_amm,
         phosphate = mean_pho)

#export to csv
#write_csv(nact_nuts1, "data/Nact/nuts/nact_nuts_full_processed.csv")
