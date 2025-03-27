# Scratch code 
summary_stats_WT <- as.data.frame(summary_stats_WT)
parameters <- rownames(summary_stats_WT[2:5,])
ggplot(summary_stats_WT, aes(x = rownames(summary_stats_WT[2:5], y = summary_stats_WT[2:5,4], fill=parameters)))

### Come back to this
# Update this section later 
# Good for visualization and summarizing all plots in one place 

# Regression data set
#####################################

rm(list=ls())
set.seed(1280)
x<-runif(50)
y<-3+2*x+rnorm(50, mean=0, sd=2)

dat<-data.frame(y, x)


# Bayesian regression in R+WinBUGS
# or R+OpenBUGS
#####################################


library(arm)
library(rjags)
library(coda)

N<-length(x)

dat<- list("N" = N, "y" = y, "x" = x)
jags.inits <- function() {list (alpha = 0, beta = 0, tau.y = 1)}
parameters<- c("alpha", "beta", "tau.y")

reg.jags<-jags.model(file="example1.bug", data=dat, inits=jags.inits, n.chains=1, n.adapt=1000)
update(reg.jags, n.iter=1000) # burn in
regression.sim<-coda.samples(reg.jags, variable.names=parameters, n.iter=15000)

geweke.diag(regression.sim)
heidel.diag(as.mcmc(regression.sim))
raftery.diag(as.mcmc(regression.sim))

summary(regression.sim)

library(mcmcplots)
mcmcplot(regression.sim, dir=getwd())

denplot(regression.sim)
caterplot(regression.sim, c("alpha", "beta"), val.lim=c(-1,6))
abline(v=0., lty=2)














# Regression with Unit Heterogeneity data set
#####################################

rm(list=ls())
set.seed(129382013)

unit<-c()
x<-c()
y<-c()

for(i in 1:20){
  
  unit.t<-rep(i, 50)
  unit.effect<-rnorm(1, mean=0, sd=3)
  x.t<-runif(50)
  y.t<-3+2*x.t+rnorm(50, mean=0, sd=2)+unit.effect
  
  unit<-c(unit, unit.t)
  x<-c(x, x.t)
  y<-c(y, y.t)
  
}

dat<-data.frame(y, x, unit)



# Bayesian RE regression in R+WinBUGS
# or R+OpenBUGS
#####################################

library(arm)
library(rjags)
library(coda)


N<-length(x)
M<-length(unique(unit))

dat<- list("N" = N, "M" = M, "unit" = unit, "y" = y, "x" = x)
jags.inits <- function() {list (alpha = 0, beta = 0, tau.y = 1, tau.re = 0.5, ranef.v = rep(0, 20))}
parameters<- c("alpha", "beta", "tau.y", "tau.re")


re.mod.jags<-jags.model(file="example3.bug", data=dat, inits=jags.inits, n.chains=4, n.adapt=1000)
update(re.mod.jags, n.iter=1000) # burn in
re.sim<-coda.samples(re.mod.jags, variable.names=parameters, n.iter=20000, n.chain=4)


geweke.diag(re.sim)
heidel.diag(re.sim)
raftery.diag(re.sim)

library(mcmcplots)
mcmcplot(re.sim, dir=getwd())

denplot(re.sim)
denplot(re.sim, collapse=T)
caterplot(re.sim, c("alpha", "beta"), val.lim=c(-1,6))
abline(v=0, lty=2)








library(dclone)
library(snow)
load.module("lecuyer")

cl<-makeCluster(10, type="SOCK")
parLoadModule(cl, "lecuyer", quiet=T)

dat<- list("N" = N, "M" = M, "unit" = unit, "y" = y, "x" = x)
jags.inits <- function() {list (alpha = 0, beta = 0, tau.y = 1, tau.re = 0.5, ranef.v = rep(0, 20))}
parameters<- c("alpha", "beta", "tau.y", "tau.re")

re.par.inits<-parallel.inits(jags.inits, n.chains=10)
parJagsModel(cl, name="par.re.mod", file="example3.bug", data=dat, inits=re.par.inits, n.chains=10, n.adapt=1000)

parUpdate(cl, "par.re.mod", n.iter=1000) # burn-in
par.re.sim<-parCodaSamples(cl, "par.re.mod", variable.names=parameters, n.iter=8000, n.chain=10)

head(par.re.sim)

denplot(par.re.sim)
denplot(par.re.sim, collapse=T)

stopCluster(cl)




##################################################################################
# Come back to this: quality check model 
# Load necessary libraries
library(arm)
library(rjags)
library(coda)
library(mcmcplots)
load("~/Sinai/Sweis Lab/Projects/Shank3/Outputs/Run1/M1_params_WT - Copy.RData")

# Step 3: MCMC Diagnostics - Check convergence and autocorrelation
# Geweke diagnostic: Assesses convergence by comparing means of the beginning and end parts of the chain
geweke_diag <- geweke.diag(results_WT)
print(geweke_diag)

# Heidelberger & Welch diagnostic: Tests stationarity and convergence
heidel_diag <- heidel.diag(as.mcmc(results_WT))
print(heidel_diag)

# Raftery & Lewis diagnostic: Tests for the number of iterations required to obtain accurate results
raftery_diag <- raftery.diag(as.mcmc(results_WT))
print(raftery_diag)

# Summary of the MCMC results
summary_stats_WT2 <- summary(results_WT)
print(summary_stats_WT2)

# Step 4: Visualizations using 'mcmcplots' package
# Traceplot to visualize chains
mcmcplot(results_WT, dir = getwd())

# Density plot for each parameter
denplot(results_WT)

# Caterplot for alpha and beta (boundary separation and drift rate)
caterplot(results_WT, c("alpha.mu", "b1.mu"), val.lim = c(-1, 6))

# Add vertical lines at specific values for reference
abline(v = 0., lty = 2)

output_path = "C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3/Outputs"

# Optional: Save the visualizations to files
output_filepath <- file.path(output_path, "M1_params_WT_plots.png")
png(output_filepath)
mcmcplot(results_WT, dir = getwd())
dev.off()

# Step 5: Save results and model outputs
save(results_WT, Data_WT_attempt2, summary_stats_WT2, file = output_filepath)



##############################################################
# Come back to: Save settings 
group = "WT"
run = "R1"
device = "BM"
if (device == "FS") { 
  output_dir = path("C:/Users/Feede/Documents/Sam_Temp/Shank3/Outputs/FS_Outputs")
}  else { 
  output_dir = path("C:/Users/Sammb/Documents/Sinai/Sweis Lab/Projects/Shank3/Outputs/BM_Outputs")
}
name <- paste0(device, "_", group, "_", run)
filename <- paste0(name,".RData")
save_as <- path(paste0(output_dir,"/",filename))

# Save files 
Dataset <- Data_WT
Results <- results_WT
SummaryStats <- summary_stats_WT
save(Dataset, Results, SummaryStats, file = save_as)


# Variables to save: 
# Data_WT_attempt2
# results_WTs
# summary_stats_WT2

