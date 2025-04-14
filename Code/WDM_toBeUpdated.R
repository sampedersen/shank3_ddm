####      DDM - Shank - Batch Processing 
########################################################################################################################################################

# Reset & settings 
rm(list=ls()) 
set.seed(1) 
options(digits = 3) 

# Load packages 
pack = c("rtdists", "RWiener", "tidyr", "dplyr", "purrr", "ggplot2", "gridExtra", 
         "gt", "xml2", "reshape2", "coda", "readxl", "patchwork", "runjags", 
         "tidyverse", "DEoptim", "Rcpp", "parallel", "RcppParallel", "loo", 
         "stats4", "pracma", "tidymodels", "fdrtool", "boot", "ggridges", 
         "ggpubr", "HDInterval", "fs", "here", "readr", "writexl")

# Define paths 
setwd("C:/Users/Sammb/Documents/Addtnl Code/dmc-master/dmc-master")
source("dmc/dmc.R")

########################################################################################################################################################

# Identify the device ["BM", "FS"]
Device = "BM"     # BM or FS 
# Attempt Date
Attempt_Date = "Apr08"
# Attempt Number 
Attempt_Num = 1
# Approach ["WDM", "RJAGS"]
Approach = "WDM"  
# Data file name 
data_filename = "Epoch4_Full_DDM_Data.xlsx"  
# Identify the sub-epoch 

########################################################################################################################################################

# Set day limits 
Day_Min_Lookup <- c(`1`=18, `2`=23, `3`=28, `4`=33, `5`=38, `6`=43, `7`=48)
Day_Max_Lookup <- c(`1`=22, `2`=27, `3`=32, `4`=37, `5`=42, `6`=47, `7`=52)


# Set home dir and source location based on device 
if (Device == "BM"){
  setwd("C:/Users/Sammb/Documents/Addtnl Code/dmc-master/dmc-master")
  source("dmc/dmc.R")
  home_dir = path("C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3")
  dedicated_cores <- round((detectCores(logical=TRUE)-2.5)/2)
} else {
  home_dir = path("C:/Users/Feede/Documents/Sam_Temp/Shank3")
  setwd("C:/Users/Feede/Documents/Sam_Temp/Shank3/Code/dmc-master/dmc-master")
  source("dmc/dmc.R")
  dedicated_cores <- round((detectCores(logical=TRUE)-2.5)/2)}
RcppParallel::setThreadOptions(numThreads = dedicated_cores)

# Set output location depending on device and model approach 
if (Approach == "WDM"){
  output_path = home_dir / "Outputs/WDM"
} else {
  output_path = home_dir / "Outputs/RJAGS"
}

# Set up paths 
data_path = home_dir / "Data"
data.pathname = home_dir / "Data"
model_path = home_dir / "Model"
code_path = home_dir / "Code"
session_name = paste(Approach,Device,Attempt_Date,Attempt_Num,sep="_")

########################################################################################################################################################

# Identify the full, overall dataset 
shank_file <- file.path(data.pathname, data_filename)
mice_data_raw <- read_excel(shank_file)
mice_data_raw <- as.data.frame(mice_data_raw)
numMice <- length(unique(mice_data_raw$mouse))             

# Select for variables of interest  
mice_data = mice_data_raw %>%
  rename(
    mouse = `mouse`, 
    day = `day`,
    trial = `trial`,
    choice = `OZ outcome 1 or 0`,  # Column has spaces, so must be in backticks
    rt = `offer zone RT (s)`,  
    offer = `offer`,        
    group = `genotype terminal`,
    value = `offer value`,
    rank = `pellet rank`
    ) %>%
  select(mouse, day, trial, choice, rt, offer, group, value, rank)

# Set up columns for stimulus and response type 
mice_data = mice_data %>%
  mutate(
    S = ifelse(value<0, "s1", "s2"),
    R = ifelse(choice==0,"r1","r2"),
    myType = case_when(
      S == "s1" & R == "r1" ~ "s1.r1",
      S == "s1" & R == "r2" ~ "s1.r2",
      S == "s2" & R == "r1" ~ "s2.r1",
      S == "s2" & R == "r2" ~ "s2.r2"),
    RespType = case_when(
      myType %in% c("s1.r1", "s2.r2") ~ "1",
      myType %in% c("s1.r2", "s2.r1") ~ "0"))

# Remove RT outliers         
mice_data = mice_data %>%
  group_by(mouse) %>%
  filter(rt > max(.250,(mean(rt) - (2*sd(rt)))), 
         rt < min(7,(mean(rt) + (2*sd(rt))))) %>% 
  ungroup()

# Compute the number of trials lost due to RT filtering 
lossToCleaning <- (nrow(mice_data_raw) - nrow(mice_data)) / nrow(mice_data_raw)

# Separate for training (80%) and testing (20%)
train_test_split <- mice_data %>%            
  group_by(mouse,day) %>%
  mutate(Set = ifelse(runif(n()) < 0.8, "Train","Test"))
train_data <- filter(train_test_split,Set=="Train")
test_data <- filter(train_test_split,Set=="Test")  

# Reassign target dataset
mice_data_cleaned <- mice_data
mice_data <- train_data

# Initialize storage
wdm_storage <- list()

# Loop over mouse names
for (mouse_id in unique(mice_data$mouse)) {
  # Output which mouse is being processed 
  print(paste("Processing Mouse:", mouse_id))
  # Pull current mouse's data 
  current_mouse <- mice_data %>%
    filter(mouse == mouse_id) %>%
    mutate(
      q = rt, 
      condition = group,
      resp = ifelse(RespType == "1", "upper", "lower") %>% as.character(),
      group = factor(ifelse(S == "s1", "s1", "s2"))
    ) %>%
    select(mouse, day, q, resp, group, rank, condition)
  # Declare a list for storing data by preference 
  wdm_pref <- list()
  # Loop through rankings 
  for (j in 1:4) {
    # Filter by rank
    wdm_data <- current_mouse %>% filter(rank == j)
    
    # Filter by stim
    wdm_data_s1 <- wdm_data %>% filter(group == "s1")
    wdm_data_s2 <- wdm_data %>% filter(group == "s2")
    
    # Split each stim by epoch
    wdm_data_s1_epochs <- list()
    wdm_data_s2_epochs <- list()
    
    for (epoch in 1:7) {
      Day_Min <- Day_Min_Lookup[as.character(epoch)]
      Day_Max <- Day_Max_Lookup[as.character(epoch)]
      
      # Subset each stim group by epoch
      wdm_data_s1_epochs[[epoch]] <- wdm_data_s1 %>%
        filter(day >= Day_Min & day <= Day_Max)
      
      wdm_data_s2_epochs[[epoch]] <- wdm_data_s2 %>%
        filter(day >= Day_Min & day <= Day_Max)
    }
    
    # Store both stim groups (each with 7 epochs)
    wdm_pref[[j]] <- list(
      wdm_data_s1 = wdm_data_s1_epochs,
      wdm_data_s2 = wdm_data_s2_epochs
    )
  }
  
  # Save for this mouse
  wdm_storage[[as.numeric(mouse_id)]] <- wdm_pref
}

# Storage format : 
# wdm_storage[["1"]][[1]][["wdm_data_s1"]][[1]]
# wdm_storage[[mouse]][[rank]][["s1/s2"]][[epoch]]






# Initialize output
wdm_param_storage <- list()

# Loop over mice
for (iMice in 1:numMice) {
  param_byPref <- list()  # Reset per mouse
  
  for (iRank in 1:4) {
    params_byStim <- list()
    
    for (stim_name in c("wdm_data_s1", "wdm_data_s2")) {
      params_byEpoch <- list()
      
      for (iEpoch in 1:7) {
        wdm_data <- wdm_storage[[(iMice)]][[iRank]][[stim_name]][[iEpoch]]
        
        # Fit only if non-empty
        if (nrow(wdm_data) > 0) {
          fit <- wdm(wdm_data)  # Your model function
        } else {
          fit <- NA
        }
        
        # Store with epoch index
        params_byEpoch[[paste0("epoch", iEpoch)]] <- fit
      }
      
      # Store stimulus type
      params_byStim[[stim_name]] <- params_byEpoch
    }
    
    # Store by rank
    param_byPref[[paste0("rank", iRank)]] <- params_byStim
  }
  
  # Store by mouse
  wdm_param_storage[[as.character(iMice)]] <- param_byPref
}


long_format_params <- bind_rows(
  lapply(seq_along(wdm_param_storage), function(i) {
    bind_rows(
      lapply(seq_along(wdm_param_storage[[i]]), function(j) {
        bind_rows(
          lapply(seq_along(wdm_param_storage[[i]][[j]]), function(s) {
            bind_rows(
              lapply(seq_along(wdm_param_storage[[i]][[j]][[s]]), function(e) {
                
                # Extract coefficients for each model fit
                coeffs <- wdm_param_storage[[i]][[j]][[s]][[e]]$coefficients
                
                # Correctly access the condition from wdm_storage
                condition <- wdm_storage[[i]][[j]][[s]][[e]]$condition
                
                # Create a data frame where condition is repeated for each coefficient
                data.frame(
                  Mouse = rep(i, length(coeffs)),
                  Preference = rep(j, length(coeffs)),
                  Stimulus = rep(s, length(coeffs)),
                  Epoch = rep(e, length(coeffs)),
                  Parameter = names(coeffs),  # Coefficient names
                  Value = coeffs,  # Coefficient values
                  Condition = rep(condition, length(coeffs)),  # Repeat condition for each coefficient
                  row.names = NULL  # Avoid row names warning
                )}))}))}))}))

# Remove row names
rownames(long_format_params) <- NULL


# Export
file_name <- "parameters_est_6.csv"
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








