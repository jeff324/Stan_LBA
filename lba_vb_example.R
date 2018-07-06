rm(list=ls())
setwd('C:/Users/Jeffrey/Dropbox/backup/Research Projects/Rstan/Code/supplementary_code/LBA/')
source('lba-math.r')
require('rstan')
require('msm')

#simulation parameters
TEST_LENGTH = 200
NUM_SUBJ = 20
NUM_COND = 3
NUM_CHOICES = 2

#mcmc parameters
model = 'lba_hier.stan'


#group-level parameters
A_mu = .5
k_mu = .5
t0_mu = .5
v1_mu_c1 = 2
v2_mu_c1 = 2
v3_mu_c1 = 2
v1_mu_c2 = 3
v2_mu_c2 = 2
v3_mu_c2 = 1
v1_mu_c3 = 4
v2_mu_c3 = 2
v3_mu_c3 = 3
s = 1

#individual-level parameters
A = rtnorm(NUM_SUBJ,A_mu,.5,0,Inf)
k = rtnorm(NUM_SUBJ,k_mu,.5,0,Inf)
t0 = rtnorm(NUM_SUBJ,t0_mu,.5,0,1)
v1_c1 = rtnorm(NUM_SUBJ,v1_mu_c1,1,0,Inf)
v2_c1 = rtnorm(NUM_SUBJ,v2_mu_c1,1,0,Inf)
v3_c1 = rtnorm(NUM_SUBJ,v3_mu_c1,1,0,Inf)
v1_c2 = rtnorm(NUM_SUBJ,v1_mu_c2,1,0,Inf)
v2_c2 = rtnorm(NUM_SUBJ,v2_mu_c2,1,0,Inf)
v3_c2 = rtnorm(NUM_SUBJ,v3_mu_c2,1,0,Inf)
v1_c3 = rtnorm(NUM_SUBJ,v1_mu_c3,1,0,Inf)
v2_c3 = rtnorm(NUM_SUBJ,v2_mu_c3,1,0,Inf)
v3_c3 = rtnorm(NUM_SUBJ,v3_mu_c3,1,0,Inf)

#make simulated data 
RT = array(NA,dim = c(NUM_SUBJ,NUM_COND,2,TEST_LENGTH))
for(i in 1:NUM_SUBJ){
     out1 = rlba(TEST_LENGTH,k[i]+A[i],A[i],c(v1_c1[i],v2_c1[i]),1,t0[i])
     out2 = rlba(TEST_LENGTH,k[i]+A[i],A[i],c(v1_c2[i],v2_c2[i]),1,t0[i])
     out3 = rlba(TEST_LENGTH,k[i]+A[i],A[i],c(v1_c3[i],v2_c3[i]),1,t0[i])
     RT[i,1,1,] = out1$rt
     RT[i,2,1,] = out2$rt
     RT[i,3,1,] = out3$rt
     RT[i,1,2,] = out1$resp
     RT[i,2,2,] = out2$resp
     RT[i,3,2,] = out3$resp
}

#make data list to give to Stan
data = list(RT=RT, TEST_LENGTH=TEST_LENGTH, NUM_SUBJ=NUM_SUBJ, NUM_COND=NUM_COND, NUM_CHOICES=NUM_CHOICES)

#compile the stan model
model_cmp = stan_model(file = model,verbose = T)

#run the model using ADVI
fit <- vb(object = model_cmp,iter=50000)

#model summary
print(fit)

