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

