# Header  ======================================================================
# Title: WDM_wip.R
# Description: 
#       
# Contributors: 

# Date Created: Apr 3 2025

# Dependencies: 
#       - R packages listed in line [XX]; use line [xx] to install and line [XX]
#             to load as needed 

# Usage: [ WIP ]
#
#_______________________________________________________________________________
#=========================       Notes     =====================================
#_______________________________________________________________________________

#_______________________________________________________________________________
#====================== 1. Set up environment  =================================
#_______________________________________________________________________________
# Clear objects from environment 
rm(list = ls())
# Set random seed for reproducibility
b_seed = 10051990
set.seed(b_seed)
# Install and load packages 
pack = c("rtdists", "RWiener", "tidyr", "dplyr", "purrr", "ggplot2", "gridExtra", 
         "gt", "xml2", "reshape2", "coda", "readxl", "patchwork", "runjags", 
         "tidyverse", "DEoptim", "Rcpp", "parallel", "RcppParallel", "loo", 
         "stats4", "pracma", "tidymodels", "fdrtool", "boot", "ggridges", 
         "ggpubr", "HDInterval", "fs", "here", "readr", "writexl")
# Optional - install packages if not already installed 
new_packages <- pack[!(pack %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# Load required packaged (suppress start-up messages) 
suppressPackageStartupMessages(lapply(pack, require, character.only = TRUE))

#_______________________________________________________________________________
#====================== 1a. User-Defined Settings ==============================
#_______________________________________________________________________________
# Identify the device 
Device = "BM"     # BM or FS 
# Condition 
Condition = "HE"  # HE or WT 
# Attempt Date
Attempt_Date = "Apr03"
# Attempt Number 
Attempt_Num = 1
# Approach
Approach = "WDM"  # WDM or RJAGS

# Data file name 
data_filename = "Epoch4_Mid-to-Late_DDM_Data.xlsx"     

#_______________________________________________________________________________
#=============== 2. Set up directory paths & import data =======================
#_______________________________________________________________________________
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

# Set up paths 
data_path = home_dir / "Data"
model_path = home_dir / "Model"
output_path = home_dir / "Outputs"
code_path = home_dir / "Code"
session_name = paste(Approach,Device,Condition,Attempt_Date,Attempt_Num,sep="_")

# Load data
# Note: edit this section depending on how the data is stored 
data_filepath = file.path(data_path, data_filename)    
df_raw = read_excel(data_filepath)                    # use this for .xlsx

#_______________________________________________________________________________
#==========================    3. Pre-process data     =========================
#_______________________________________________________________________________
# Pull needed variables, rename for conventions 
df = df_raw %>%
  mutate(
    mouse = `mouse`,  # Ensure column name matches exactly
    day = `day`,
    trial = `trial`,
    choice = `OZ outcome 1 or 0`,  # Column has spaces, so must be in backticks
    rt = `offer zone RT (s)`,  
    offer = `offer`,        
    group = `genotype terminal`,
    value = `offer value`,
    S = ifelse(value < 0, "s1", "s2"),
    R = ifelse(choice == 0, "r1", "r2"),
    myType = case_when(
      S == "s1" & R == "r1" ~ "s1.r1",
      S == "s1" & R == "r2" ~ "s1.r2",
      S == "s2" & R == "r1" ~ "s2.r1",
      S == "s2" & R == "s2" ~ "s2.r2"),
    RespType = case_when(
      myType == "s1.r1" | myType == "s2.r2" ~ "1",
      myType == "s1.r2" | myType == "s2.r1" ~ "0")) %>%
  select(mouse, day, trial, choice, rt, offer, group, value, S, R, myType, RespType)

# Remove RT Outliers         
# Note: Hard limit cut-offs at [.25s,7s]
# Note: Relative limits +/- 2SD of average RT per subject 
df = df %>%
  group_by(mouse) %>%
  filter(rt > max(.250,(mean(rt) - (2*sd(rt)))), 
         rt < min(7,(mean(rt) + (2*sd(rt))))) %>% 
  ungroup()

# Compute the number of trials lost due to RT filtering 
# Note: ~5% loss in data is fine, shouldnt exceed ~8% 
lossToCleaning <- (nrow(df_raw) - nrow(df)) / nrow(df_raw)


# Separate for training (80%) and testing (20%)
train_test_split <- df %>%            # Create a new dataframe 
  group_by(mouse,day) %>%          # Group by mouse and by day 
  mutate(Set = ifelse(runif(n()) < 0.8, "Train","Test"))
train_data <- filter(train_test_split,Set=="Train")   # Assign to training dataframe
test_data <- filter(train_test_split,Set=="Test")     # Assign to testing dataframe 

# Reassign df as training_data for development and troubleshooting purposes 
# Change this later, if needed 
df_training <- train_data

# Storage for DDM info 
wdm_storage <- list()
wdm_param_storage <- list()





