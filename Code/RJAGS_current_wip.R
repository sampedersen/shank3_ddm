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
Condition = "WT"  # HE or WT 
# Attempt Date
Attempt_Date = "Apr04"
# Attempt Number 
Attempt_Num = 1

# Model Settings 
BurnIn = 1000
Sample = 10000
Thinning = 10
Chains = 2




#_______________________________________________________________________________
#=============== 2. Set up directory paths & import data =======================
#_______________________________________________________________________________

if (Device == "BM"){
  home_dir = path("C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3")
  num_cores <- detectCores(logical=TRUE)
  dedicated_cores <- round((num_cores-1.5)/2)
} else {
  home_dir = path("C:/Users/Feede/Documents/Sam_Temp/Shank3")
  num_cores <- detectCores(logical=TRUE)
  dedicated_cores <- round((num_cores-1.5)/2)
}

# Limit core usage for parallel processing 
RcppParallel::setThreadOptions(numThreads = dedicated_cores)

# Set up paths 
data_path = home_dir / "Data"
model_path = home_dir / "Model"
output_path = home_dir / "Outputs"
code_path = home_dir / "Code"

# Load data
# Note: edit this section depending on how the data is stored 
data_filename = "Epoch4_Mid-to-Late_DDM_Data.xlsx"     
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
         group = `genotype terminal`) %>%    # variable participant's condition ('WT'=Wildtype, 'HT'=Shank3 het)
  select(subj_idx, group, day, trial, choice, rt, offer)

# Summarize data prior to processing (optional) 
# Plot the number of trials per subject 
df_summary <- df %>%
  group_by(subj_idx) %>%
  summarise(n = length(rt)) %>%
  as.data.frame()
ggplot(df_summary, aes(x = n)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Number of Trials per Subject (Pre Cleaning)",x = "Number of Trials",y = "Frequency") +
  theme_minimal()
# RT distributions
RT_extremes <- c(min(df$rt), max(df$rt))
ggplot(df,aes(x=rt))+
  geom_histogram(binwidth = .5, fill = "skyblue", color = "black", alpha = 0.7)+
  labs(title="RT Distributions (Pre-Cleaning)", x = "Rts",y = "Freq") +
  theme_minimal() 

# Remove RT Outliers         
# Note: Hard limit cut-offs at [.25s,7s]
# Note: Relative limits +/- 2SD of average RT per subject 
df = df %>%
  group_by(subj_idx) %>%
  filter(rt > max(.250,(mean(rt) - (2*sd(rt)))), 
         rt < min(7,(mean(rt) + (2*sd(rt))))) %>% 
  ungroup()

# Summarize data after processing (optional) 
# Plot the number of trials per subject 
df_summary <- df %>%
  group_by(subj_idx) %>%
  summarise(n = length(rt)) %>%
  as.data.frame()
ggplot(df_summary, aes(x = n)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Number of Trials per Subject (Post Cleaning)",x = "Number of Trials",y = "Frequency") +
  theme_minimal()
# RT distributions
RT_extremes <- c(min(df$rt), max(df$rt))
ggplot(df,aes(x=rt))+
  geom_histogram(binwidth = .05, fill = "skyblue", color = "black", alpha = 0.7)+
  labs(title="RT Distribs (Post Cleaning)",x = "Rts",y = "Freq") +
  theme_minimal() 

# Compute the number of trials lost due to RT filtering 
# Note: ~5% loss in data is fine, shouldnt exceed ~8% 
lossToCleaning <- (nrow(df_raw) - nrow(df)) / nrow(df_raw)

# Make rejected choices into negative RT values 
idx <- which(df$choice == 0)  # Pull indices for rejected offer trials 
df$RT <- df$rt                # Duplicate RTs to a new column   
df$RT[idx] <- df$rt[idx] * -1 # Update new column for signed RT values 

# Separate for training and testing 
# Note: 80% of trials per day per mouse should be used for training model,
#       remaining 20% should be saved for testing model performance 
train_test_split <- df %>%            # Create a new dataframe 
  group_by(subj_idx,day) %>%          # Group by mouse and by day 
  mutate(Set = ifelse(runif(n()) < 0.8, "Train","Test"))
# Notes for function above: 
#     - Create a new column, Set, using mutate() function
#     - runif() generates random numbers from a uniform distribution between 0 and 1 
#     - n() is the number of rows of data (number of trials for each mouse, each day)
#     - runif(n()) generates a number of random values based on n()
#     - If the value generated is less than 0.8, assign to Training set; otherwise, 
#         assign to Testing set 
#     - Should result in approximately 80/20 train/test split 
train_data <- filter(train_test_split,Set=="Train")   # Assign to training dataframe
test_data <- filter(train_test_split,Set=="Test")     # Assign to testing dataframe 

# (Optional) Write training/testing datasets out as .xlsx
training_filename = "training_dataset.xlsx"
training_filepath = data_path / training_filename
testing_filename = "testing_dataset.xlsx"
testing_filepath = data_path / testing_filename
write_xlsx(train_data, training_filepath)
write_xlsx(test_data, testing_filepath)

### TEMPORARY: Reassign df as training_data for development and troubleshooting purposes 
# Change this later, if needed 
df <- train_data

#_______________________________________________________________________________
#==================== # 4. Model fitting (DEVELOPMENT PHASE) ===================
#_______________________________________________________________________________
# 4. Model fitting (DEVELOPMENT PHASE) 

# Note: Fitting model for WT group solely right now
# Note: Code can be adapted to for-loop over each group, but can be computationally 
#       expensive
# Note: For now, do WT and evaluate before attempting HT group 
# Note: Uses Model 1 (no drift intercept)

# 4a. Prepare data for WT group ------------------------------------------------

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
ns=length(unique(idxP))     # number of subjects 

# Prepare data in JAGS format 
dat <- dump.format(list(N=N, 
                        y=y, 
                        idxP=idxP, 
                        offer=offer, 
                        rt=rtpos, 
                        ns=ns))

#   4b. Parameter initialization -----------------------------------------------
# Notes on variables: 
#     - Initializing parameters for MCMC Chains 
#     - Each parameter has a mean value (.mu) - average across subjects
#     - And precision parameter (.pr or bias.kappa) - how much individual subjects vary around the mean 
#         - Higher value - more precision/less individual variability around mean 
#         - Lower value - less precision/more individual variability around mean 
#     - (TEMP) Able to initialize 3 sets, but for DEVELOPMENT PHASE, only using 1 

# Variable meanings and assignments 
#      - alpha.mu - mean boundary separation (decision threshold)
#      - alpha.pr - precision for boundary separation 

#      - theta.mu - mean non-decision time (Ter)
#      - theta.pr - precision for non-decision time

#      - b1.mu - mean drift rate effect of offer
#      - b1.pr - precision for drift rate effect 

#      - bias.mu - mean starting point bias (initial preference)
#         - bias.mu = 0.5  ==> no bias
#         - bias.mu < 0.5  ==> bias towards rejecting
#         - bias.mu > 0.5  ==> bias towards accepting 
#      - bias.kappa - precision for bias

#      - y_pred - initial predicted RT values 

#      - .RNG.name - specifies random number generator 
#         - inits3 - Super-Duper
#         - inits2 - Wichmann-Hill
#         - inits1 - Mersenne-Twister
#       - .RNG.seed - random seed (99999) for reproducibility 

# Specify model parameters to monitor 
#     - Deviance = how well the model fits the data (minimize as much as possible)
#     - Group level parameters - overall effects to boundary, Ter, bias, drift effects 
#     - Subject level parameters - individual differences in alpha, theta, bias, b1
monitor = c("deviance",
            # Group-level parameters (boundary, non-decision time, starting-point bias, effect of offer on drift rate)
            "alpha.mu","theta.mu","bias.mu","b1.mu",
            # Subject-level parameters
            "alpha.p","theta.p","bias.p","b1.p")

# Run the JAGS model using M1_ug_drift.txt and parallel MCMC sampling 
# - Specify 3 MCMC chains (how many to run in parallel; increasing is more 
#     robust but slower, while lower is faster but may converge poorly 
# - Plots posterior distributions automatically (plots=true)
# Robust / speed tradeoffs 
#     - Burnin: number of samples to discard; removes transient effects 
#           - Increasing is more stable/slower, decreasing is faster/poorer convergence
#     - Sample: Number of iterations to retain for analysis 
#           - Increasing is more stable/slower, decreasing is faster/poorer convergence
#     - Thin: Sampling frequency for reducing auto-correlations
#           - Increasing is faster/poorer convergence, decreasing is more stable/slower 
#           - Keeping every 20th iteration, vs every 10th iteration 
# Initial parameters 1:
inits1 <- dump.format(list(
  alpha.mu=2, alpha.pr=.5,            # Alpha
  theta.mu=0.200, theta.pr=0.005,     # Theta
  b1.mu=0.16, b1.pr=0.1,              # Drift rate
  bias.mu=0.52, bias.kappa=.5,        # Bias
  y_pred=y,                           # RT values
  .RNG.name="base::Super-Duper",      # Random Number Generator
  .RNG.seed=67882))                   # Seed selected 

# Initial parameters 2: 
inits2 <- dump.format(list(
  alpha.mu=1, alpha.pr=0.5,         # Alpha 
  theta.mu=0.2,theta.pr=0.005,         # Theta 
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
Results <- run.jags(model = file.path(model_path,"M1_ug_drift.txt"), 
                    monitor=monitor, data=dat, n.chains=Chains,
                    inits=c(inits1,inits2, inits3), 
                    plots = TRUE, method="parallel", module="wiener",
                    burnin=BurnIn,
                    sample=Sample,
                    thin=Thinning)

# Save summary statistics
Summary<-summary(Results)

# Save model outputs and results for later analysis
# Generates file name based on group
output_filename = paste(Device,Condition,Attempt_Date,Attempt_Num, sep="_")
output_filename = "M1_params_HE.RData"
output_filepath = output_path / output_filename
save(results_WT2,Data_WT_attempt2,summary_stats_WT2, file=output_filepath) 


#_______________________________________________________________________________
#=====================    # 4c. Sample Console Outputs      ====================
#_______________________________________________________________________________
# "Calling the simulation using the parallel method..."
# "Following the progress of chain 1 (the program will wait for all chains to finish before continuing):"
# Model is being ran in parallel using multiple CPU cores 
# Console is showing progress for only one chain, but JAGS is running both chains
#     simultaneously 
# You'll only see progress for one chain at a time, but all chains must finish 
#     before proceeding 

# "Welcome to JAGS 4.3.1 on Thu Mar 20 14:43:31 2025"
# "JAGS is free software and comes with ABSOLUTELY NO WARRANTY"
# Confirms which version of JAGS is running with timestamp for when model fitting began 

# "Loading module: basemod: ok"
# Core JAGS module 
# "Loading module: bugs: ok"
# Required module for compatibility with WinBUGS-style models 
# ". Loading module: wiener: ok"
# Required for DDM 
# ". Loading module: dic: ok"
# Required for Deviance Information Criterion (DIC) calculations for model comparison


#". . Reading data file data.txt"
# JAGS is reading the data file (data.txt), which was created internally by run.jags
#" . Compiling model graph"
# JAGS is translating the model into a computational graph 
#"    Resolving undeclared variables"
# JAGS checks that all variables in the model are defined 
#"    Allocating nodes"
# JAGS sets up probability distributions and dependencies 

# "Graph information:"
#"    Observed stochastic nodes: 85113"
# Number of datapoints 
#" Unobserved stochastic nodes: 170241"
# Unknown parameters that JAGS is estimating 
#"    Total graph size: 520885"
# Total number of elements in the model

#" . Reading parameter file inits1.txt"
# Loading in the initial parameter values from inits1.txt
# Even though inits3 is specified, run.jags often generates multiple initsX.txt files 
#" . Initializing model"
# Initialize with starting values 

#" . Adapting 1000"
# -------------------------------------------------| 1000
# ++++++++++++++++++++++++++++++++++++++++++++++++++ 100%
# Adaptation successful
# JAGS adapts the MCMC sampler for 1000 iterations 
# Successful adaptation means JAGS tuned its samplers for efficient sampling 

# . Updating 2000
# -------------------------------------------------| 2000
# **************************** 
# Burn-in phase begins, with 2000 iterations 
# Model is running without storing results to allow the chains to converge 
# JAGS will begin collecting posterior samples after 

#  . . . . . . . . . . Updating 10000
#  -------------------------------------------------| 10000
# Sampling phase 
# JAGS is running 10,000 iterations across all chains 
# Only keeps every 20 of the iterations 

# ==============================================================================
# 5. Check diagnostics 
#-------------------------------------------------------------------------------  
# Note: This point and beyond is still a WIP and pending model development and outputs 




################################################################################
# Check convergence of MCMC chains 
# Check using the Potential Scale Reduction Factor (PSRF) 
# - PSRF < 1.1 generally acceptable indication of chain convergence 
# - PSRF > 1.1 generally indicates lack of proper convergence and requires either
#     more sample or burn-in 
# Note: These were hard coded; need to double check and see what the outputs are
#   each saved as and add proper naming 


# Load and check WT group 
gg = "WT"   
load(file = output_filepath)
summary_stats_WT2 %>% as.data.frame() %>%     # Evaluate chain convergence 
  filter(psrf > 1.1)              # Value should be <1.1 for all params

# Load and check HE group (***DOUBLE CHECK NAMING CONVENTIONS)
gg = "HE"
load(file = output_path / paste0("M1_Params_",gg,".RData"))
summary_stats_HE2 %>% as.data.frame() %>%
  filter(psrf > 1.1)

################################################################################
# Extract, analyze, & compare parameters

# Initialize params in an empty dataframe 
params = NULL

# Loop through the two group, WT and HE 
for (gg in c("WT", "HE")) {
  # Merge MCMC chains
  # For each group, load the MCMC chain results 
  load(file = output_path / paste0("M1_Params_",gg,".RData"))
  # Combine the MCMC chains from the three different chains stored within the results 
  # rbind() stacks the results from the three chains into a single dataframe/matrix 
  chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
  # Pull unique subject IDs 
  subj = unique(Data$idxP)
  
  # Loop through each subject 
  for (s in subj){
    # For each subject, extract mouse's data from Data object 
    subj_idx = unique(Data$subj_idx[Data$idxP == s])
    # Extract the parameters for each subject from the chain 
    # Parameters include regression weight for drift rate, boundary separation, 
    #   non-decision time, starting bias 
    # Calculate the parameters as an average across the MCMC chains 
    # Store parameters in a temporary dataframe (tmp_res) 
    tmp_res = data.frame(
      idxP = s,
      subj_idx = subj_idx,
      group = gg,
      wOffer = mean(chain[,c( paste( c("b1.p[",toString(s),"]"), collapse = ""))]),
      boundary = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))]),
      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))]),
      bias = mean(chain[,c( paste( c("bias.p[",toString(s),"]"), collapse = ""))])
    )
    params = rbind(params, tmp_res)
  }
}

# Within Groups t-tests (Pairwise Comparisons)
params %>%
  # Reshape the dataframe to long format 
  # Transforms the columns into a single column (name) and corresponding value 
  #   in adjacent column (value)
  pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
  # Group the data by their group (WT or HE) and the parameter name (wOffer, boundary, etc)
  group_by(group,name) %>%
  # Compute a paired t-test to compare parameter values between two sessions,
  #   store the tvalue and pvalue of t-test 
  summarise(tvalue = t.test(value ~ session, paired = TRUE)$statistic,
            pvalue = t.test(value ~ session, paired = TRUE)$p.value) %>%
  # Only show parameters where p-value is less than 0.05 (indicates significant 
  #     differences between sessions for those parameters)
  filter(pvalue < .05)  

# Between Groups t-test (Between-Group Comparison)
# Compare parameter values between groups (WT and HE) across sessions, rather than 
#     within a session 
params %>%
  # Make into a long format
  pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
  # Mutate to turn the groups into factors with two levels (WT, HE)
  mutate(group = factor(group, levels = c("WT","HE"))) %>%
  # Group by session and parameter name 
  group_by(session,name) %>%
  # Compute t-test to compare values of each parameter between groups for each session
  summarise(tvalue = t.test(value ~ group)$statistic,
            pvalue = t.test(value ~ group)$p.value) %>%
  # Only show parameters where p-value is less than 0.05 (indicates significant 
  #     differences between groups for each session for those parameters)
  filter(pvalue < .05)

# One-way ANOVA for between-group comparison (Post-Hoc Tukey Test)
# Perform one-way ANOVA on wOffer parameter to test is there are differences between
#   the two groups 
summary(aov(data = params,
            formula = wOffer ~ group))
# If the ANOVA shows significant differences between groups, Tukey Honestly 
#   Significant Difference (HSD) post-hoc test identifies which specific groups are 
#   significantly different from each other 
TukeyHSD(aov(data = params,
             formula = wOffer ~ group),which = "group")

################################################################################
# Visualization 
ggplot(params,aes(x = group, y = wOffer, color = group, group = group)) +
  theme_pubr(base_size = 18) +
  geom_point(position = position_dodge2(width=.25), alpha = .5,size=3) +
  geom_hline(yintercept =0,linetype="dashed",color="grey",linewidth=1) +
  stat_summary(geom = "line",position = position_dodge2(width=.25),linewidth=1.5) +
  stat_summary(position = position_dodge2(width=.25),size=1.5,linewidth=1.5) +
  scale_color_brewer(type="qual",palette = 4) +
  labs(x = "Dx",
       y = "Influence of Offer on Drift Rate",
       color = "Dx")