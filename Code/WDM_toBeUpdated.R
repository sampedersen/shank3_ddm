########################################################################################################################################################
####      DDM - Shank - Batch Processing 
########################################################################################################################################################

# Workflow 
#

########################################################################################################################################################
####      Housekeeping
########################################################################################################################################################

# Reset & settings 
rm(list=ls()) 
set.seed(1) 
options(digits = 3) 

# Load packages 
library(rtdists) # DDM functions rdiffusion(), ddiffusion(), qdiffusion()
library(RWiener) # alternate DDM functions rwiener(), deviance.wdm(), logLik.wdm()
library(tidyr)   # unnest() function
library(dplyr)   # mutate() and arrange() functions
library(purrr)   # split() and map() functions
library(ggplot2)   # graphing
library(gridExtra) # for multiple ggplots in one window
library(gt)        # grammar of tables
library(xml2)      # save formatted tables as rtf 
library(reshape2)  # for melt() function
library("coda")
library(readxl)
library(patchwork)

# Define paths 
setwd("C:/Users/Sammb/Documents/Addtnl Code/dmc-master/dmc-master")
data.pathname = "C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank_DDM/Data/Processed/HDDM_setup"
source("dmc/dmc.R")

########################################################################################################################################################
####      Set up for Batching 
########################################################################################################################################################

# Identify the full, overall dataset 
output_path <- "C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3/Outputs/WDM"
shank_file <- file.path(data.pathname, "ddm_shankData_processed_longitudinal.xlsx")
miceSheets <- excel_sheets(shank_file)    # Pull mice sheet names 
numMice <- length(miceSheets)                   # number of mice based on number of sheets 
mice_data <- list()

# Import sheets to list for storage 
for (i in 1:numMice) {  
  mice_data[[i]] <- read_excel(shank_file, sheet = miceSheets[i])}

# Clean/adjust data 
for (i in 1:length(mice_data)) {  
  current_mouse <- mice_data[[i]]  
  current_mouse <- current_mouse %>%
    mutate(
      S = ifelse(Value < 0, "s1", "s2"),
      R = ifelse(Response == 1, "r1", "r2"),
      myType = case_when(
        S == "s1" & R == "r1" ~ "s1.r1",
        S == "s1" & R == "r2" ~ "s1.r2",
        S == "s2" & R == "r1" ~ "s2.r1",
        S == "s2" & R == "r2" ~ "s2.r2"),
      RespType = case_when(
        myType == "s1.r1" | myType == "s2.r2" ~ "1",
        myType == "s1.r2" | myType == "s2.r1" ~ "0")) %>%
    relocate(S, .after = Value) %>%
    relocate(R, .after = Response)
  
  # Optional: Filter RTs 
  current_mouse <- current_mouse %>%
    filter(RTs<=4)
  
  # Put back into list 
  mice_data[[i]] <- current_mouse
}

# Declare lists to store DDM info
wdm_storage <- list()
wdm_param_storage <- list()

# Generate simulated dataset
for (i in seq_along(mice_data)) {  # Iterate over mice
  current_mouse <- mice_data[[i]]  
  wdm_pref <- list()  # Reset per mouse
  for (j in 1:4) {  # Iterate over preferences
    mouse_subset <- current_mouse[current_mouse$Rankings == j,]
    # Set up training dataset
    wdm_data <- mouse_subset %>%
      dplyr::mutate(
        q = RTs, 
        resp = factor(ifelse(RespType == 1, "upper", "lower")),
        group = factor(ifelse(S == "s1", "s1", "s2"))
      ) %>%
      dplyr::select(q, resp, group) %>%
      as.data.frame()
    
    # Separate by stimulus type
    wdm_data_byStim <- list(
      wdm_data_s1 = wdm_data %>% dplyr::filter(group == "s1") %>% as.data.frame(),
      wdm_data_s2 = wdm_data %>% dplyr::filter(group == "s2") %>% as.data.frame())
    wdm_pref[[j]] <- wdm_data_byStim  # Store by preference
  }
  wdm_storage[[i]] <- wdm_pref  # Store by mouse
}

# Fit WDM model
for (i in seq_along(mice_data)) {
  param_byPref <- list()  # Reset per mouse
  for (j in 1:4) {
    params <- list()  # Reset per preference
    for (s in 1:2) {
      wdm_data <- wdm_storage[[i]][[j]][[s]]
      # Ensure the data is non-empty before fitting
      if (nrow(wdm_data) > 0) {
        params[[s]] <- wdm(wdm_data)
      } else {
        params[[s]] <- NA  # Handle empty datasets
      }
    }
    param_byPref[[j]] <- params
  }
  wdm_param_storage[[i]] <- param_byPref  # Store params by mouse
}

# Convert to long format
long_format_params <- bind_rows(
  lapply(seq_along(wdm_param_storage), function(i) {
    bind_rows(
      lapply(seq_along(wdm_param_storage[[i]]), function(j) {
        bind_rows(
          lapply(seq_along(wdm_param_storage[[i]][[j]]), function(s) {
            coeffs <- wdm_param_storage[[i]][[j]][[s]][['coefficients']]
            data.frame(
              Mouse = i,
              Preference = j,
              Stimulus = s,
              Parameter = names(coeffs),  # Store coefficient names
              Value = coeffs  # Store coefficient values
            )}))}))}))
rownames(long_format_params) <-  NULL

# Export
file_name <- "parameters_est_3.csv"
export_file <- file.path(output_path, file_name)  # Properly concatenate paths
write.csv(long_format_params, export_file, row.names = FALSE)

# Initialize lists properly
LLEs <- vector("list", length(mice_data))  
for (i in seq_along(mice_data)) {
  LLEs_pref <- vector("list", 4)  # List for preferences
  
  for (j in 1:4) {
    LLEs_stim <- vector("list", 2)  # List for stimulus types
    
    for (s in 1:2) {
      if (!is.null(wdm_param_storage[[i]][[j]][[s]])) {  
        LLEs_stim[[s]] <- wdm_param_storage[[i]][[j]][[s]][['loglik']]
      } else {
        LLEs_stim[[s]] <- NA  # Assign NA correctly
      }
    }
    
    LLEs_pref[[j]] <- LLEs_stim  # Store per preference
  }
  
  LLEs[[i]] <- LLEs_pref  # Store per mouse
}








