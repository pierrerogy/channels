# Calculate elemental composition of resources

# Load packages
library(tidyverse)

# Load data
res_raw <- 
  read.csv("data_exp2/resource_content.csv",
           sep = ";")
# Check data
str(res_raw)

# Average per resource type
res <- 
  res_raw %>% 
  dplyr::group_by(resource, unit) %>% 
  dplyr::select(-sample) %>% 
  dplyr::summarise_all(mean,
                na.rm = TRUE) 

  

# Convert mg.kg to %
res <- 
  res %>% 
  dplyr::filter(unit == '%') %>% 
  dplyr::bind_rows(res %>% 
                     dplyr::filter(unit == 'mg.kg') %>% 
                     dplyr::select(-C, -N, -S.1) %>% 
                     dplyr::mutate(across(everything(), as.numeric()),
                                   across(everything(), ~./1e+05)))

  
  dplyr::mutate(.cols
  
