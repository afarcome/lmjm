simula_joint_HMcontinuous <- function(n=100,k=3,sepxiv=1,type=c("norm","bino")){
	
# Simulation of the from joint model with continuous HM chain described in
# "A shared-parameter continuous-time hidden Markov and survival model for longitudinal data
#  with informative drop-out"
# by F.Bartolucci and A.Farcomeni (2018)
#
# INPUT:
#
# n          = sample size
# k          = number of latent states
# sepxiv     = separation between the latent state support points 
# type       = type of response variables ("norm" or "bino")
#
# OUTPUT:
#
# Y      = matrix of responses (each row corresponds to an individual)
# Y0     = matrix of all responses (without dropout)
# W      = matrix of individual covariates for the survival process
# X      = array of individual covariates on the longitudinal process
# tv     = vector of survival times
# censv  = vector of indicators of censoring
# Q      = infinitesimal transition matrix
# piv    = initial probability vector
# bev    = regression coefficients on the longitudinal process
# nu     = parameter of the Weibull distribution
# xiv    = vector of support points
# phi    = relation between longitudal and survival process
# psiv   = vector of regression paraemters on the survival process
# si2    = variane for normal responses
# TT     = time of each longitudinal observaton
# TT0    = time of each longitudinal observaton (without dropout)
# tU     = matrix of latent states (continuous time)
# tT     = matrix of continuous time jumps
# U      = matrix of latent states (discrete time)
# Q      = infinite transition matrix
# dev    = indicator for non-censoring time


# set up
	type = match.arg(type)
	cv = rep(10,n)                  # censoring times (v stands for vector)
	xiv = seq(1:k)*sepxiv           # support points
	xiv = xiv-mean(xiv)           
	piv = rep(1/k,k)                # initial probabilty vector
	Q = matrix(1/(k-1),k,k)         # infinitesimal transition matrix
	diag(Q) = -1
	if(type == "norm")	si2 = 0.25 else si2 = NULL   # variance of continuous data
	X = array(0,c(n,10+1,2))      # array of covariates (AR(1))
	for(i in 1:n){
		X[i,1,] = rnorm(2)
		for(t in 2:11) X[i,t,] = 0.9*X[i,t-1,]+rnorm(2)*sqrt(0.19)
	}
	bev = c(-1,1)                      # corresponding regression parameters 
	phi = 0.5                          # random-effect coefficient on survival process
	W = cbind(1,matrix(rnorm(2*n),n))  # covariates for survival process
	psiv = c(-4,-1,1)                  # corresponding regression parameters
	nu = 2                             # parameter of the Weibull distribution 

# simulate data
	tU = tT = tTau = matrix(NA,n,100)     # for MC process (initial t stands fro tilde)
	tjv = rep(0,n)                        # number of jumps
	R = Q; diag(R) = 0                    # jump matrix
	qtot = rowSums(R)                     # exponential sojourn rate
	R = R/qtot

# draw continuous MC process
	for(i in 1:n){
		tT[i,1] = 0
		tU[i,1] = which(rmultinom(1,1,piv)==1)
		j = 0
		while(tT[i,j+1]<=cv[i]){
			j = j+1
			tTau[i,j] = rexp(1,qtot[tU[i,1]])
			tT[i,j+1] = tT[i,j]+tTau[i,j]
			tU[i,j+1] = which(rmultinom(1,1,R[tU[i,j],])==1)
		}
		tjv[i] = j
	}

# draw survival time
	bv = tv = rep(0,n)            # vector if random draws and surivavial times
	censv = rep(FALSE,n)          # vector with indicators of censoring
	for(i in 1:n){
		Ht = rep(0,tjv[i])          # cumulative hazard function
		ind = 1:(tjv[i]-1)
		if(tjv[i]>1) Ht[2:tjv[i]] = exp(as.vector(W[i,]%*%psiv))*cumsum(exp(xiv[tU[i,ind]]*phi)*
		                            (tT[i,ind+1]^nu-tT[i,ind]^nu))
		bv[i] = -log(runif(1))
		di = sum(Ht<=bv[i])
		tv[i] = (tT[i,di]^nu + (bv[i]-Ht[di])/(exp(xiv[tU[i,di]]*phi+W[i,]%*%psiv)))^(1/nu)
		if(tv[i]>cv[i]){
			censv[i] = TRUE
			tv[i] = cv[i]
		}
	}
	dev = 1-censv                 # indicator for non-censoring time

# draw longitudinal process
	jv = rep(0,n)                           # number of observations
	TT = TT0 = Y = Y0 = matrix(NA,n,11)     # time of each observation and value
	U = matrix(NA,n,100)                    # matrix of latent state at each observation
	for(i in 1:n){
		jv[i] = floor(tv[i])+1
		TT0[i,] = 0:10
		TT[i,1:jv[i]] = 0:(jv[i]-1)
		U[i,1] = tU[i,1]
		tmp = xiv[U[i,1]]+X[i,1,]%*%bev
		if(type=="norm") Y0[i,1] = tmp+rnorm(1)*sqrt(si2)     # continuous case
		if(type=="bino"){                                     # binary case
			tmp = exp(tmp)/(1+exp(tmp))  
			Y0[i,1] = 1*(runif(1)<tmp)
		}
		for(j in 2:11){
			tj = sum(tT[i,]<=TT0[i,j],na.rm=TRUE)
			U[i,j] = tU[i,tj]
			tmp = xiv[U[i,j]]+X[i,j,]%*%bev
			if(type=="norm") Y0[i,j] = tmp+rnorm(1)*sqrt(si2)
			if(type=="bino"){
				tmp = exp(tmp)/(1+exp(tmp))
				Y0[i,j] = 1*(runif(1)<tmp)
			}
		}
		Y[i,1:jv[i]] = Y0[i,1:jv[i]]
	}

# final output
return(list(Y=Y,Y0=Y0,W=W,X=X,tv=tv,censv=censv,Q=Q,piv=piv,bev=bev,nu=nu,xiv=xiv,
	          phi=phi,psiv=psiv,si2=si2,TT=TT,TT0=TT0,tU=tU,tT=tT,U=U,Q=Q,dev=dev))

}
