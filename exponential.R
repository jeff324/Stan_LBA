rm(list=ls())
setwd("C:/...") #needs to be same directory as exponential.stan
library(rstan)
#make simualated data
dat = rexp(500,1)
len = length(dat)
#run the Stan model
fit <- stan(file   = "exponential.stan", 
            data   = list(Y=dat,LENGTH=len),          
            warmup = 750,                 
            iter   = 1500,                
            chains = 3)                   
#model summary
print(fit)
#posterior predictions
mcmc_chain = as.matrix(fit)
lambda = mean(mcmc_chain[,'lambda'])
pred = mcmc_chain[,'pred']
hist(dat,probability=T)
lines(density(pred))
#autocorrelation plot
acf(mcmc_chain[,'lambda'])
#samples vs. iteration plot
traceplot(fit)
#plot lambda
layout(mat=matrix(1:2,1,2))
par(mar=c(4,4,.1,.5))
hist(dat,probability=T,main='',xlab='y')
box()
lines(density(mcmc_chain[,'pred']),lwd=2)
plot(density(mcmc_chain[,'lambda']),main='',xlab=expression(lambda),las=1)
