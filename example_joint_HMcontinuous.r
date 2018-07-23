# Example of application of the R functions related to the paper:
#
# "A shared-parameter continuous-time hidden Markov and survival model for longitudinal data
#  with informative drop-out"
# by F.Bartolucci and A.Farcomeni (2018)

# Simulation of data (see the function code for the help)

# Normal distribution for the responses
source("simula_joint_HMcontinuous.R")
out = simula_joint_HMcontinuous()

# Estimate with complete data
tT = out$tT   # time of any (time-continuous jump)
TT = out$TT   # observation times
tU = out$tU   # latent states
TT = out$TT   # time of observation of any longitudinal outcome
W = out$W     # matrix of baseline covariates
X = out$X     # array of time-varying covariates
Y = out$Y     # array of time-varying covariates
tv = out$tv   # vector of survival times
dev = out$dev # vector of censoring
source("est_comp_joint_HMcontinuous.R")
est = est_comp_joint_HMcontinuous(tT,tU,TT,W,X,Y,tv,dev,k=3,type="norm")

# Estimate with incomplete data
source("est_incomp_joint_HMcontinuous.R")
est1 = est_incomp_joint_HMcontinuous(Y,TT,X,W,tv,dev,k=3,type="norm",verbose=TRUE,stderr=TRUE)

# Bernoulli distribution for the responses
out1 = simula_joint_HMcontinuous(n=500,type="bino")

# Estimate under with complete data
tT = out1$tT   # time of any (time-continuous jump)
TT = out1$TT   # observation times
tU = out1$tU   # latent states
TT = out1$TT   # time of observation of any longitudinal outcome
W = out1$W     # matrix of baseline covariates
X = out1$X     # array of time-varying covariates
Y = out1$Y     # array of time-varying covariates
tv = out1$tv   # vector of survival times
dev = out1$dev # vector of censoring
est2 = est_comp_joint_HMcontinuous(tT,tU,TT,W,X,Y,tv,dev,k=3,type="bino")

# Estimate with incomplete data
source("est_incomp_joint_HMcontinuous.R")
est3 = est_incomp_joint_HMcontinuous(Y,TT,X,W,tv,dev,k=3,type="bino",verbose=TRUE,stderr=TRUE)
