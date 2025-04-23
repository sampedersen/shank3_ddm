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
b_seed = 100
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
#====================== 1a. Model settings ====== ==============================
#_______________________________________________________________________________
BurnIn = 10000 # First 10k samples we toss
Sample = 2000  # Sample 2000 times
Thinning = 10  # 1 per every 10 samples picked 
Chains = 6    # 3 starting points (min) to converge at same end point 
# number of cores to recruit (1:1 core:chains)
RcppParallel::setThreadOptions(Chains)

# Data file name 
data_filename = "Epoch4_Mid-to-Late_DDM_Data.xlsx"     

#_______________________________________________________________________________
#====================== 2. Directory Set-Up ===================================
#_______________________________________________________________________________
# Specify the device depending on user account directories that exist 
if (dir.exists("C:/Users/Sammb/")){
  Device = "BM"
  home_dir = path("C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3")
  setwd(home_dir)
} else if(dir.exists("C:/Users/Feede")){
  Device = "FS"
  home_dir = path("C:/Users/Feede/Documents/Sam_Temp/Shank3")
  setwd(home_dir)
} else {print("Device directory not recognized.")}

# Set up paths 
data_path = home_dir / "Data"
model_path = home_dir / "Model"
output_path = home_dir / "Outputs"
code_path = home_dir / "Code"

# Load data
data_filepath = file.path(data_path, data_filename)    
df_raw = read_excel(data_filepath)                    


# Attempt Date
Attempt_Date = format(Sys.Date(), "%B-%d")
# Attempt Number 
Attempt_Num = format(Sys.time(), "%H-%M")


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
         group = `genotype terminal`, # variable participant's condition ('WT'=Wildtype, 'HT'=Shank3 het)
         value = `offer value`) %>%    
  select(subj_idx, group, day, trial, choice, rt, offer, value)

# Remove RT Outliers (note: removed bottom limit, all RTs are under 2SDs)         
df = df %>%
  group_by(subj_idx) %>%
  filter(rt < (mean(rt) + (2*sd(rt)))
         #rt < (mean(rt) + (2*sd(rt)))
         ) %>% 
  ungroup()

# Compute the number of trials lost due to RT filtering 
lossToCleaning <- (nrow(df_raw) - nrow(df)) / nrow(df_raw)
print(paste0("Loss to cleaning: ", lossToCleaning*100, "%"))

# Make rejected choices into negative RT values 
idx <- which(df$choice == 0)  # Pull indices for rejected offer trials 
df$RT <- df$rt                # Duplicate RTs to a new column   
df$RT[idx] <- df$rt[idx] * -1 # Update new column for signed RT values 

#_______________________________________________________________________________
#==================== # 4. Model fitting (DEVELOPMENT PHASE) ===================
#_______________________________________________________________________________

# Subsetting 
set.seed(22)
selected_mice <- df %>%
  distinct(subj_idx,group) %>%    # Lists each mouse with condition info 
  group_by(group) %>%   # Group by condition 
  slice_sample(n=12) %>%  # 12 per group 
  pull(subj_idx)    # Pull subject ids 
print(sort(selected_mice))
selected_days <- c(23, 28, 32, 37, 42, 47, 52)
subset_df <- df %>%
  filter(subj_idx %in% selected_mice) %>%
  filter(day %in% selected_days)

subset_df <- subset_df %>%
  mutate(value_group = if_else(value<0,"neg","pos"))
balanced_subset <- subset_df %>%
  group_by(subj_idx,choice,value_group) %>%
  sample_frac(.25) %>%
  ungroup()

nrow(balanced_subset %>% filter(choice==0))
nrow(balanced_subset %>% filter(choice==1))
nrow(balanced_subset %>% filter(value_group=="pos"))
nrow(balanced_subset %>% filter(value_group=="neg"))











# 4. Model fitting (DEVELOPMENT PHASE) 
for (genotype in unique(balanced_subset$group)) {
  Data = balanced_subset %>%
    filter(group==genotype) %>%
    mutate(idxP = as.numeric(ordered(subj_idx)))
    
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
  
  monitor = c("deviance",            # Come back to this 
              "alpha.mu","theta.mu","bias.mu","b1.mu", # Group-level parameters (boundary, non-decision time, starting-point bias, effect of offer on drift rate)
              "alpha.p","theta.p","bias.p","b1.p")    # Subject-level parameters
  
  inits1 <- dump.format(list(
    alpha.mu=1.9, 
    alpha.pr=.6,            # Alpha
    theta.mu=0.05, 
    theta.pr=0.5,       # Theta
    b1.mu=0.30, 
    b1.pr=0.5,              # Drift rate
    bias.mu=0.75, 
    bias.kappa=.25,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Super-Duper",      # Random Number Generator (Blair has had issues w this one (not so super))
    .RNG.seed=10002
    ))                            # Seed selected 
  
  inits2 <- dump.format(list(
    alpha.mu=1.5, 
    alpha.pr=.25,            # Alpha
    theta.mu=0.150, 
    theta.pr=0.3,       # Theta
    b1.mu=0.01, 
    b1.pr=0.20,              # Drift rate
    bias.mu=0.5, 
    bias.kappa=.5,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Wichmann-Hill", 
    .RNG.seed=67882
  ))  
  
  inits3 <- dump.format(list(
    alpha.mu=2.6, 
    alpha.pr=.05,            # Alpha
    theta.mu=0.14, 
    theta.pr=0.7,       # Theta
    b1.mu=0.7, 
    b1.pr=0.3,              # Drift rate
    bias.mu=0.25, 
    bias.kappa=.75,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Mersenne-Twister", 
    .RNG.seed=6666))
  
  inits4 <- dump.format(list(
    alpha.mu=2.00, 
    alpha.pr=0.50,            # Alpha
    theta.mu=0.130, 
    theta.pr=0.5,       # Theta
    b1.mu=0.16, 
    b1.pr=0.10,              # Drift rate
    bias.mu=0.52, 
    bias.kappa=.3,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Super-Duper",      # Random Number Generator (Blair has had issues w this one (not so super))
    .RNG.seed=999
  ))                            # Seed selected 
  
  inits5 <- dump.format(list(
    alpha.mu=1.0, 
    alpha.pr=.25,            # Alpha
    theta.mu=0.120, 
    theta.pr=0.25,       # Theta
    b1.mu=0.05, 
    b1.pr=0.40,              # Drift rate
    bias.mu=0.5, 
    bias.kappa=.65,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Wichmann-Hill", 
    .RNG.seed=333
  ))  
  
  inits6 <- dump.format(list(
    alpha.mu=2.5, 
    alpha.pr=.75,            # Alpha
    theta.mu=0.110, 
    theta.pr=0.75,       # Theta
    b1.mu=0.1, 
    b1.pr=0.6,              # Drift rate
    bias.mu=0.60, 
    bias.kappa=.45,        # Bias
    y_pred=y,                           # RT values
    .RNG.name="base::Mersenne-Twister", 
    .RNG.seed=222))
  


  # Store results 
  Results <- run.jags(model = file.path(model_path,"M1_ug_drift.txt"), 
                      monitor=monitor, data=dat, n.chains=Chains,
                      inits=c(inits1,inits2, inits3, inits4,inits5,inits6), 
                      plots = TRUE,
                      method="parallel", 
                      modules ="wiener",  
                      burnin=BurnIn,
                      sample=Sample,
                      thin=Thinning)
  
  # Save summary statistics
  Summary<-summary(Results)
  
  # Save model outputs and results for later analysis
  output_filename = paste(paste(genotype,Attempt_Date,Attempt_Num, sep="_"), ".RData",sep="")
  output_filepath = output_path / output_filename
  save(Results,Data,Summary, file=output_filepath) 
}

