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
Chains = 3    # 3 starting points (min) to converge at same end point 
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

# Attempt Date
Attempt_Date = "April-23"
# Attempt Number 
Attempt_Num = "13-43"


# Check diagnostics

# Load one of the files
for (gg in c("WT", "HE")) {
  if (gg=="WT") {
    load(file = paste0(output_path,"/WT_",Attempt_Date,"_",Attempt_Num,".RData"))
    WT_Results <- Results
    WT_Summary <- Summary
    WT_Data <- Data
    rm(Results, Summary, Data)
    filtered <- WT_Summary %>% as.data.frame() %>% filter(psrf > 1.1)
    print(paste("WT - problematic psrf entries:", nrow(filtered)))
    print(filtered)
  } else if (gg=="HE") {
    load(file=paste0(output_path,"/HE_",Attempt_Date,"_",Attempt_Num,".RData"))
    HE_Results <- Results
    HE_Summary <- Summary
    HE_Data <- Data
    rm(Results, Summary, Data)
    filtered <- HE_Summary %>% as.data.frame() %>% filter(psrf > 1.1)
    print(paste("HE - problematic psrf entries:", nrow(filtered))) 
    ## If there are any > 1.1, will need to re-run with more burn-in
    print(filtered)
  }}

# Extract parameters
params = NULL
for (gg in c("WT", "HE")) {
  if (gg=="WT") {
    load(file = "C:/Users/Feede/Documents/Sam_Temp/Shank3/Outputs/April-24/WT_April-24_09-26.RData")
    for (i in 1:length(Results[["mcmc"]])) {
      # If chain is empty, initialize it with the first element
      if (is.null(chain)) {
        chain <-Results$mcmc[[i]]
      } else {
        # Otherwise, bind the new element to the chain
        chain <- rbind(chain, Results$mcmc[[i]])
  } else if (gg=="HE") {
    load(file="C:/Users/Feede/Documents/Sam_Temp/Shank3/Outputs/April-24/WT_April-24_09-26.RData")
    for (i in 1:length(HE_Results[["mcmc"]])) {
      # If chain is empty, initialize it with the first element
      if (is.null(chain)) {
        chain <- WT_Results$mcmc[[i]]
      } else {
        # Otherwise, bind the new element to the chain
        chain <- rbind(chain, WT_Results$mcmc[[i]])
      }
  
    }
  
  
  
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
param_filename <- file.path(output_path, "params.xlsx")
write_xlsx(params,param_filename)

params %>%
  # Reshape data to long format
  # Take columns for wOffer, boundary, nDT, and bias
  # Convert into two columns (name and value) 
  pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
  group_by(group,name) %>%
  summarise(tvalue = t.test(value ~ session, paired = TRUE)$statistic,
            pvalue = t.test(value ~ session, paired = TRUE)$p.value) %>%
  filter(pvalue < .05)

params %>%
  pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
  mutate(group = factor(group, levels = c("HC","BD"))) %>%
  group_by(group,name) %>%
  summarise(tvalue = t.test(value ~ group)$statistic,
            pvalue = t.test(value ~ group)$p.value) %>%
  filter(pvalue < .05)

summary(aov(data = params,
            formula = wOffer ~ group))
TukeyHSD(aov(data = params,
             formula = wOffer ~ group),which = "group")


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