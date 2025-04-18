# Header  ======================================================================
# Title: RJAGS_Shank3_DDM.R
# Description: 
#       Fit a Drift Diffusion Model (DDM) to model choice behaviors in offer zone
#       of Restaurant Row for wildtype (WT) and Shank3 haploidinsufficient (HT) 
#       mice.
# Contributors: Blair Shevlin (primary code author), Sam Pedersen (modified and
#     adapted for current use)
# Date Created: Mar 20 2025
# Dependencies: 
#       - R packages listed in line [XX]; use line [xx] to install and line [XX]
#             to load as needed 
# Usage: [ WIP ]
#
#_______________________________________________________________________________
#=========================       Notes     =====================================
#_______________________________________________________________________________
# RJAGS is an R package for utilizing the Just Another Gibbs Sampler (JAGS) software
# JAGS is a program used for performing Bayesian statistical modeling via 
# Markov Chain Monte Carlo (MCMC) methods 

#_______________________________________________________________________________
#====================== 1. Set up environment  =================================
#_______________________________________________________________________________
# Clear objects from environment 
rm(list = ls())
# Set random seed for reproducibility
b_seed = 10051990
set.seed(b_seed)

# Install +/- load packages 
pack = c("runjags","tidyverse","purrr","DEoptim","Rcpp","parallel","RcppParallel","loo","coda",
         "stats4","pracma","tidymodels","fdrtool","boot","ggridges","ggpubr",
         "HDInterval","fs","here","readxl","dplyr","readr","writexl")
# Optional - install packages if not already installed 
new_packages <- pack[!(pack %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# Load required packaged (suppress start-up messages) 
suppressPackageStartupMessages(lapply(pack, require, character.only = TRUE))

#' Define helper function for standard error (se)
#' @param x A numeric vector for which the standard error is to be computed.
#' @return A numeric value representing the standard error of the input vector.
#' @examples
#' # Example usage:
#' se(c(1, 2, 3, 4, 5))
se <- function(x) {
  sd(x) / sqrt(length(x))}

#_______________________________________________________________________________
#====================== 1a. User-Defined Settings ==============================
#_______________________________________________________________________________
# Identify the device 
Device = "FS"     # BM or FS 
# Condition 
Condition = "HE"  # HE or WT 
# Attempt Date
Attempt_Date = format(Sys.Date(), "%B-%d")
# Attempt Number 
Attempt_Num = format(Sys.time(), "%H-%M")

# Model Settings 
# Red flags lol

# Randomly sample a subset of the dataset instead of reducing sampling/thinning
# Blair ~95k burnin
# Increase 
BurnIn = 10000 # First 10k samples we toss
Sample = 2000  # Sample 2000 times
Thinning = 10  # 1 per every 10 samples picked 
Chains = 3    # 3 starting points (min) to converge at same end point 

# Data file name 
data_filename = "Epoch4_Mid-to-Late_DDM_Data.xlsx"     

#_______________________________________________________________________________
#=============== 2. Set up directory paths & import data =======================
#_______________________________________________________________________________
# Limit number of cores used for parallel processing based on number available per device 
if (Device == "BM"){
  home_dir = path("C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3")
  setwd(home_dir)
} else {
  home_dir = path("C:/Users/Feede/Documents/Sam_Temp/Shank3")
  setwd(home_dir)
}
RcppParallel::setThreadOptions(Chains)

# Set up paths 
data_path = home_dir / "Data"
model_path = home_dir / "Model"
output_path = home_dir / "Outputs"
code_path = home_dir / "Code"

# Load data
# Note: edit this section depending on how the data is stored 
# Check if file exists as .csv, if not check as .xlsx
data_filepath = file.path(data_path, data_filename)    
df_raw = read_excel(data_filepath)                    # use this for .xlsx

#_______________________________________________________________________________
#==========================    3. Pre-process data     =========================
#_______________________________________________________________________________
# Pull needed variables, rename for conventions 
df = df_raw %>%
  mutate(subj_idx = `mouse`,     # variable containing mouse number 
         day = `day`,            # variable containing day of testing 
         trial = `trial`,        # variable containing trial numbers
         choice = `OZ outcome 1 or 0`,       # variable containing choice (coded, 1=accept, 0=reject)
         rt = `offer zone RT (s)`,           # variable containing rt in seconds
         offer = `offer`,        # variable containing offer amount
         group = `genotype terminal`,
         value = `offer value`) %>%    # variable participant's condition ('WT'=Wildtype, 'HT'=Shank3 het)
  select(subj_idx, group, day, trial, choice, rt, offer, value)

# Remove RT Outliers         
# Note: Hard limit cut-offs at [.25s,7s]
# Note: Relative limits +/- 2SD of average RT per subject 

# Consideration: Bin the response time acrose subjects (15,25ms) and see the timepoint that 
# they're "correct" around 50%
# Easy case where mice knows its preference, stronger true decision
# Floor can def be lower than .25, read mice DDMs 

###### COME BACK 
df = df %>%
  group_by(subj_idx) %>%
  filter(rt > (mean(rt) - (2*sd(rt))), 
         rt < (mean(rt) + (2*sd(rt)))) %>% 
  ungroup()

# Compute the number of trials lost due to RT filtering 
# Note: ~5% loss in data is fine, shouldnt exceed ~8% 
lossToCleaning <- (nrow(df_raw) - nrow(df)) / nrow(df_raw)
print(paste0("Loss to cleaning: ", lossToCleaning*100, "%"))

# Make rejected choices into negative RT values 
idx <- which(df$choice == 0)  # Pull indices for rejected offer trials 
df$RT <- df$rt                # Duplicate RTs to a new column   
df$RT[idx] <- df$rt[idx] * -1 # Update new column for signed RT values 

# TAKE THIS OUT - DONT NEED RN 
# Separate for training (80%) and testing (20%)
#train_test_split <- df %>%            # Create a new dataframe 
#  group_by(subj_idx,day) %>%          # Group by mouse and by day 
#  mutate(Set = ifelse(runif(n()) < 0.8, "Train","Test"))
#train_data <- filter(train_test_split,Set=="Train")   # Assign to training dataframe
#test_data <- filter(train_test_split,Set=="Test")     # Assign to testing dataframe 

# Reassign df as training_data for development and troubleshooting purposes 
# Change this later, if needed 
#df <- train_data

#_______________________________________________________________________________
#==================== # 4. Model fitting (DEVELOPMENT PHASE) ===================
#_______________________________________________________________________________
# Subsetting 
selected_mice <- df %>%
  distinct(subj_idx,group) %>%    # Lists each mouse with condition info 
  group_by(group) %>%   # Group by condition 
  slice_sample(n=5) %>%  # 5 per group 
  pull(subj_idx)    # Pull subject ids 
subset_df <- df %>%
  filter(subj_idx %in% selected_mice)
samples_per_trialCombo <- 50
balanced_subset <- subset_df %>%
  group_by(subj_idx,choice,value) %>%
  filter(n() >= samples_per_trialCombo) %>% # Keep groups that have enough data 
  slice_sample(n=samples_per_trialCombo)
  














# 4. Model fitting (DEVELOPMENT PHASE) 

# Filter based on user-identified group of interest 
if(Condition == "HE"){
  Data = df %>%
    filter(group=="HE") %>%
    mutate(idxP=as.numeric(ordered(subj_idx)))
} else {
  Data = df %>%
    filter(group=="WT") %>%
    mutate(idxP=as.numeric(ordered(subj_idx)))
}

Data = Data
  idxP = Data$idxP         # Sequentially numbered list of subject indices 
  offer = Data$offer       # Pull offers 
  rtpos = Data$rt          # Pull RTs (original, non-signed)
  
  # Modeling variables
  y=Data$RT                # Signed RT values (to be predicted) 
  N= length(y)                # Number of total trials
  ns=length(unique(idxP))     # number of subjects (loops through for subject-level)
  
  # Prepare data in JAGS format 
  dat <- dump.format(list(N=N, 
                          y=y, 
                          idxP=idxP, 
                          offer=offer, 
                          rt=rtpos, 
                          ns=ns))
  
#   4b. Parameter initialization -----------------------------------------------

  # Identify variables to track 
  monitor = c("deviance",            # Come back to this 
              # Group-level parameters (boundary, non-decision time, starting-point bias, effect of offer on drift rate)
              "alpha.mu","theta.mu","bias.mu","b1.mu",
              # Subject-level parameters
              "alpha.p","theta.p","bias.p","b1.p")

  # Run the JAGS model using M1_ug_drift.txt and parallel MCMC sampling 
  # Tweak these last as needed, address model "settings" first 
  # Change seed > values > RNG
  # Keep the precisions > 0.01 (ideally ~ish 0.1~0.5)
  
  # Uniform priors and uninformed priors 
  # Initial parameters 1:
  
  inits1 <- dump.format(list(
    alpha.mu=2, alpha.pr=.5,            # Alpha
    theta.mu=0.200, theta.pr=0.5,       # Theta
    b1.mu=0.16, b1.pr=0.1,              # Drift rate
    bias.mu=0.52, bias.kappa=.5,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Super-Duper",      # Random Number Generator (Blair has had issues w this one (not so super))
    .RNG.seed=67882))                   # Seed selected 
  
  # Initial parameters 2: 
  inits2 <- dump.format(list(
    alpha.mu=1, alpha.pr=0.5,           # Alpha 
    theta.mu=0.2,theta.pr=0.5,          # Theta 
    b1.mu=0.05, b1.pr=0.1,              # Drift rate 
    bias.mu=0.5,bias.kappa=.5, 
    y_pred=y,  
    .RNG.name="base::Wichmann-Hill", 
    .RNG.seed=1234))
  
  # Initial parameters 3: 
  inits3 <- dump.format(list(
    alpha.mu=2.5, alpha.pr=0.5, 
    theta.mu=0.3, theta.pr=0.05, 
    b1.mu=0.1, b1.pr=0.01, 
    bias.mu=0.6,bias.kappa=1, 
    y_pred=y, 
    .RNG.name="base::Mersenne-Twister", 
    .RNG.seed=6666))
  
  # Store results 
  Results <- run.jags(model = file.path(model_path,"M1_ug_drift.txt"), 
                      monitor=monitor, data=dat, n.chains=Chains,
                      inits=c(inits1,inits2, inits3), 
                      plots = TRUE, method="parallel", module="wiener",  ## Homework: Read about WIener process and 2d space
                      burnin=BurnIn,
                      sample=Sample,
                      thin=Thinning)

  # Save summary statistics
  Summary<-summary(Results)
  
  # Goal for 
  
  # Save model outputs and results for later analysis
  # Generates file name based on group
  output_filename = paste(Device,Condition,Attempt_Date,Attempt_Num, sep="_")
  output_filename = "M1_params_HE.RData"
  output_filepath = output_path / output_filename
  save(results_WT2,Data_WT_attempt2,summary_stats_WT2, file=output_filepath) 
