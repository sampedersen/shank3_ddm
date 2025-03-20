# ==============================================================================
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
# ==============================================================================
#           Notes
#-------------------------------------------------------------------------------
# RJAGS is an R package for utilizing the Just Another Gibbs Sampler (JAGS) software
# JAGS is a program used for performing Bayesian statistical modeling via 
# Markov Chain Monte Carlo (MCMC) methods 



# Clear objects from environment 
rm(list = ls())

# Packages to load; suppressed messages 
#install.packages("writexl")
pack = c("runjags","tidyverse","purrr","DEoptim","Rcpp","parallel","RcppParallel","loo","coda",
         "stats4","pracma","tidymodels","fdrtool","boot","ggridges","ggpubr",
         "HDInterval","fs","here","readxl","dplyr","readr","writexl")
suppressPackageStartupMessages(lapply(pack, require, character.only = TRUE))


# **** Notes needed ****


# Limits multi-threading to 1 core for parallel computing 
RcppParallel::setThreadOptions(numThreads = 1)

# Helper function for computing standard error of a numeric vector 
se = function(x) sd(x) / sqrt(length(x))

# Random seed for reproducibility
b_seed = 10051990
set.seed(b_seed)

# Define paths for data and model files 
home_dir = path("C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank_DDM")
data_path = home_dir / "Data"
preproc_data_path = data_path / "Pre-Processed"
model_path = home_dir / "Scripts/Model"
output_path = home_dir / "Outputs"


################################################################################

# Load data
#df_raw = read.csv(data_path / "Pre-Processed/Epoch4_Mid-to-Late_DDM_Data.csv")
df_raw = read_excel(preproc_data_path / "Epoch4_Mid-to-Late_DDM_Data.xlsx")

# !! Note this assume participant is a receiver in the UG !!

# Pre-processing - identify the column names associated with each variable 
df = df_raw %>%
  mutate(subj_idx = `mouse`,     # variable containing participant IDs
         day = `day`,            # variable containing day 
         trial = `trial`,        # variable containing trial numbers
         choice = `OZ outcome 1 or 0`,       # variable containing choice (coded, 1=accept, 0=reject)
         rt = `offer zone RT (s)`,           # variable containing rt in seconds
         offer = `offer`,        # variable containing offer amount
         group = `genotype terminal`) %>%    # variable participant's condition (BD=binge-drinker, HC=healthy control)
  select(subj_idx, group, day, trial, choice, rt, offer)

# Summarize data prior to processing 

# Number of trials per subject 
df_summary <- df %>%
  group_by(subj_idx) %>%
  summarise(n = length(rt)) %>%
  as.data.frame()
# plot 
ggplot(df_summary, aes(x = n)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Number of Trials per Subject (Pre Cleaning)",
       x = "Number of Trials",
       y = "Frequency") +
  theme_minimal()

# RT distributions
min(df$rt)
max(df$rt)
ggplot(df,aes(x=rt))+
  geom_histogram(binwidth = .5, fill = "skyblue", color = "black", alpha = 0.7)+
         labs(title="RT Distribs (Pre Cleaning)",
              x = "Rts",
              y = "Freq") +
         theme_minimal()

# Remove outliers 
df = df %>%
  # Hard limits at [.25s,7s]
  # Anything outside of 2 SD of mean for that subject 
  group_by(subj_idx) %>%
  filter(rt > max(.250,(mean(rt) - (2*sd(rt)))), # hard floor of 250 ms
         rt < min(7,(mean(rt) + (2*sd(rt))))) %>% # hard ceiling of 7 sec
  ungroup()

# Summarize data prior to processing 

# Number of trials per subject 
df_summary <- df %>%
  group_by(subj_idx) %>%
  summarise(n = length(rt)) %>%
  as.data.frame()
# plot 
ggplot(df_summary, aes(x = n)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Number of Trials per Subject (Post Cleaning)",
       x = "Number of Trials",
       y = "Frequency") +
  theme_minimal()

# RT distributions
min(df$rt)
max(df$rt)
ggplot(df,aes(x=rt))+
  geom_histogram(binwidth = .05, fill = "skyblue", color = "black", alpha = 0.7)+
  labs(title="RT Distribs (Post Cleaning)",
       x = "Rts",
       y = "Freq") +
  theme_minimal()

# Compute post-filtering data loss (should be ~5%, no more than 8%)
(nrow(df_raw) - nrow(df)) / nrow(df_raw)


# Separate 80% for training, 20% for testing 
# Pull 80% of trials per day per mouse 
train_test_split <- df %>% 
  group_by(subj_idx,day) %>%
  # Create a new column, Set using mutate() function
  # runif() generates random numbers from a uniform distribution between 0 and 1 
  # n() is the number of rows of data (number of trials for each mouse, each day)
  # runif(n()) generates a number of random values based on n()
  # runif(n()) generates n() random numbers from uniform distribution 
  # If the value generated is less than 0.8, assign to Training set; otherwise, 
  #   assign to Testing set 
  # Should result in approximately 80% of trials (per mouse per day / all trials)
  #   being split into training (20% testing) 
  mutate(Set = ifelse(runif(n()) < 0.8, "Train","Test"))
train_data <- filter(train_test_split,Set=="Train")
test_data <- filter(train_test_split,Set=="Test")

# Count the number of "Train" and "Test" for each subject and day
train_test_count <- train_test_split %>%
  group_by(subj_idx, day, Set) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the proportion of "Train" vs "Test"
ggplot(train_test_count, aes(x = subj_idx, y = count, fill = Set)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  facet_wrap(~day, scales = "free_x") +  # Facet by day
  labs(title = "Train vs Test Distribution per Subject per Day",
       x = "Subject",
       y = "Number of Trials",
       fill = "Set") +
  scale_fill_manual(values = c("Train" = "skyblue", "Test" = "lightcoral")) +  # Custom colors for Train and Test
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability

# Write out 
write_xlsx(train_data, preproc_data_path / "train_data.xlsx")
write_xlsx(test_data, preproc_data_path / "test_data.xlsx")

# Reassign df as training dataset for now 
df <- train_data

################################################################################
# Model fitting (Model 1: No drift intercept)
# Fit the model separately for each group 
# For each group (BD or HC)
#for (gg in unique(df$group)) {
    # Filter Data for one group at a time 
  gg = "WT"
    Data_WT = df %>%
        filter(group == gg) %>%
        mutate(idxP = as.numeric(ordered(subj_idx)))  # Convert participant IDs into numeric indices 
    # Differentiate RTs for choice type 
    idx = which(Data$choice==0)       # Identify indices rejected trials 
    Data$RT <- Data$rt                # Pull RTs to new column 
    Data$RT[idx] = Data$rt[idx] * -1  # Make rejected choices' RTs negative values 

    
    idxP = Data$idxP          # Make a sequentially numbered subj index
    offer = Data$offer        # Pull offers 
    rtpos = Data$rt           # Pull RTs (original positives)

    # Modeling variables 
    y= Data$RT                # Signed RTs
    N = length(y)             # Number of total trials 
    ns = length(unique(idxP)) # Number of subjects 

    #--------------------------------------------
    
    # Prepare data in JAGS format for processing 
    dat <- dump.format(list(N=N, y=y, idxP=idxP, offer=offer, rt=rtpos, ns=ns))

    # Initialize 3 parameters for Markov Chain Monte Carlo (MCMC) chains 
    # Each parameter has a mean value (.mu) - average across subjects
    # precision parameter (.pr or bias.kappa) - how much individual subjects vary around the mean 
        #     Higher value - more precision/less individual variability around mean 
        #     Lower value - less precision/more individual variability around mean 
    # alpha.mu - mean boundary separation (decision threshold)
    # alpha.pr - precision for boundary separation 
    
    # theta.mu - mean non-decision time (Ter)
    # theta.pr - precision for non-decision time
    
    # b1.mu - mean drift rate effect of offer
    # b1.pr - precision for drift rate effect 
    
    # bias.mu - mean starting point bias (initial preference)
    #     bias.mu = 0.5  ==> no bias
    #     bias.mu < 0.5  ==> bias towards rejecting
    #     bias.mu > 0.5  ==> bias towards accepting 
    # bias.kappa - precision for bias
    
    # y_pred - initial predicted RT values 
    
    # .RNG.name - specifies random number generator 
    #       inits3 - Super-Duper
    #       inits2 - Wichmann-Hill
    #       inits1 - Mersenne-Twister
    # .RNG.seed - random seed (99999) for reproducibility 
    inits3 <- dump.format(list( alpha.mu=1, alpha.pr=0.5, 
                                theta.mu=0.200, theta.pr=0.25,  
                                b1.mu=0.25, b1.pr=0.20, 
                                bias.mu=0.5, bias.kappa=1, 
                                y_pred=y,  
                                .RNG.name="base::Super-Duper", .RNG.seed=99999))


    inits2 <- dump.format(list( alpha.mu=2.5, alpha.pr=0.1, theta.mu=0.2,
                                theta.pr=0.05, b1.mu=0.0, b1.pr=0.01, bias.mu=0.5,
                                bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))

    inits1 <- dump.format(list( alpha.mu=2.0, alpha.pr=0.1, theta.mu=0.3,
                                theta.pr=0.05, b1.mu=-0.1, b1.pr=0.01, bias.mu=0.6,
                                bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
    
    # Specify model parameters to monitor 
    # Deviance = how well the model fits the data (minimize as much as possible)
    # Group level parameters - overall effects to boundary, Ter, bias, drift effects 
    # Subject level parameters - individual differences in alpha, theta, bias, b1 
    monitor = c("deviance",
                # Group-level parameters (boundary, non-decision time, starting-point bias, effect of offer on drift rate)
                "alpha.mu","theta.mu","bias.mu","b1.mu",
                # Subject-level parameters
                "alpha.p","theta.p","bias.p","b1.p")

    # Run the JAGS model with parallel MCMC sampling 
    # Run the JAGS model using M1_ug_drift.txt
    # Specify 3 MCMC chains 
    #   Number to run in parallel; increasing is more robust but slower, while lower
    #     is faster but may converge poorly 
    # Plots posterior distributions automatically (plots=true)
    results_WT <- run.jags(model=file.path(model_path,"M1_ug_drift.txt"), 
                        monitor=monitor, data=dat, n.chains=3, 
                        inits=inits3,
                        #inits=c(inits1,inits2, inits3), 
                        plots = TRUE, method="parallel", module="wiener", 
                        # Number of samples to burn-in/discard to remove transient effects (60k)
                        # Increasing is more stable but longer, decreasing is faster
                        # but may have poor convergence 
                        #burnin=60000, 
                        burnin=5000, 
                        # Number of iterations retained for analysis 
                        # Higher is more robust, lower is faster but noisier 
                        sample=1000, 
                        #sample=10000, 
                        # Only keep every 10th sample to reduce autocorrelations 
                        thin=10)
    # Summary statistics 
    WT_suuum<-summary(results)
    
    # Save model outputs and results for later analysis 
    # Generates file name based on group 
    save(results_WT,Data_WT,suuum_WT, file= output_path / paste0("M1_Params_",gg,".RData")) 
    
}

################################################################################
# Check diagnostics

# Load and check BD group 
gg = "BD";
load(file = output_path /  paste0("M1_Params_",gg,".RData"))
suuum %>% as.data.frame() %>%     # Evaluate chain convergence 
  filter(psrf > 1.1)              # Value should be <1.1 for all params

# Load and check HC group 
gg = "HC"
load(file = output_path / paste0("M1_Params_",gg,".RData"))
suuum %>% as.data.frame() %>%
  filter(psrf > 1.1)

## Note: if any parameters have a convergence greater than 1.1, 
##    rerun with more burn-in

################################################################################
# Extract & analyze parameters
params = NULL

for (gg in c("BD", "HC")) {
    # Merge MCMC chains
    load(file = output_path / paste0("M1_Params_",gg,".RData"))
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    subj = unique(Data$idxP)
    
    for (s in subj){
      subj_idx = unique(Data$subj_idx[Data$idxP == s])
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

# Perform t-tests to compare parameters across groups 
params %>%
  pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
  group_by(group,name) %>%
  summarise(tvalue = t.test(value ~ session, paired = TRUE)$statistic,
            pvalue = t.test(value ~ session, paired = TRUE)$p.value) %>%
  filter(pvalue < .05)
params %>%
  pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
  mutate(group = factor(group, levels = c("HC","BD"))) %>%
  group_by(session,name) %>%
  summarise(tvalue = t.test(value ~ group)$statistic,
            pvalue = t.test(value ~ group)$p.value) %>%
  filter(pvalue < .05)

summary(aov(data = params,
            formula = wOffer ~ group))
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