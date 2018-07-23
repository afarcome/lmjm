est_comp_joint_HMcontinuous <- function(tT,tU,TT,W,X,Y,tv,dev,k,type=c("norm","bino")){

# Estimate the joint model with continuous HM chain described in
# "A shared-parameter continuous-time hidden Markov and survival model for longitudinal data
#  with informative drop-out"
# by F.Bartolucci and A.Farcomeni (2018)
#
# COMPLETE DATA CASE
#
# INPUT:
# tT   = time of any (time-continuous jump)
# tU   = matrix of latent states
# TT   = matrix of times of observation of any longitudinal outcome
# W    = matrix of baseline covariates
# X    = array of time-varying covariates
# tv   = vector of survival times
# dev  = vector of censoring
# k    = number of latent states
#
# OUTPUT:
# bev  = estimated regression coefficients on the longitudinal process
# nu   = estimated parameter of the Weibull distribution
# xiv  = estimated vector of support points
# phi  = estimated relation between longitudal and survival process
# psiv = estimated vector of regression paraemters on the survival process
# piv  = estimated initial probability vector
# Q    = estimated infinitesimal transition matrix
# si2  = estimated variane for normal responses
# tnv  = vector of frequencies of each state at the first occasion
# tN   = vector of frequencies for each pair of states
# sv   = vector of sojourn times

# Preliminaries
	type = match.arg(type)
	n = nrow(tT)                    # sample size
	tjv = rowSums(!is.na(tU))-1     # number of jumps
	ncovW = ncol(W)
	TTX = dim(X)[2]
	ncovX = dim(X)[3]
# Initial values
	xiv = seq(1:k)                  # support points
	xiv = xiv-mean(xiv)           
	piv = rep(1/k,k)                # initial probabilty vector
	Q = matrix(1/(k-1),k,k)         # infinitesimal transition matrix
	diag(Q) = -1
	if(type == "norm")	si2 = 1 else si2 = NULL   # variance of continuous data
	bev = rep(0,ncovX)              # corresponding regression parameters 
	phi = 0                         # random-effect coefficient on survival process
	psiv = rep(0,ncovW)             # corresponding regression parameters
	nu = 2                          # parameter of the Weibull distribution 
	eta = log(nu)                   # reparametrization
# Compute first log-likelihood component
	tnv = rep(0,k)                  # frequency of each state at the first occasion
	for(u in 1:k) tnv[u] = sum(tU[,1]==u)
	tN = matrix(0,k,k)              # frequencies for each pair of states
	for(i in 1:n) for(j in 1:tjv[i]) tN[tU[i,j],tU[i,j+1]] = tN[tU[i,j],tU[i,j+1]]+1
	sv = rep(0,k)                   # sojourn times
	for(i in 1:n) for(j in 1:tjv[i]) sv[tU[i,j]] = sv[tU[i,j]]+tT[i,j+1]-tT[i,j]
	qv = -diag(Q)
	lk1 = sum(tnv*log(piv))
	for(u in 1:k) lk1 = lk1+sum(tN[u,-u]*log(Q[u,-u]))
	lk1 = lk1-sum(sv*qv)
# Compute second and third log-likelihood components
	par = c(bev,eta,xiv,phi,psiv)
	lk23 = lk_comp_joint_HM(par,si2,tT,tU,TT,W,X,Y,tv,dev,k,type)
# Display initial log-likelihood
	lk = lk1+lk23
	print("Initial log-likelihood")
	print(lk)  
# Estimate parameters
	piv = tnv/n
	Q = (1/sv)*tN
	qv = rowSums(Q)
	est = optim(par,lk_comp_joint_HM,gr=dlk_comp_joint_HM,method="BFGS",si2=si2,tT=tT,tU=tU,
                TT=TT,W=W,X=X,Y=Y,tv=tv,dev=dev,k=k,type=type,
                control=list(fnscale=-1,trace=1))
    par = est$par
	bev = par[1:ncovX]; eta = par[ncovX+1]; nu = exp(eta)
	xiv = par[(ncovX+2):(k+ncovX+1)]; phi = par[k+ncovX+2]; psiv = par[(k+ncovX+3):length(par)]
# estimate si2
	if(type=="norm"){
		jv = rowSums(!is.na(Y))
		si2 = 0
		for(i in 1:n){
			for(j in 1:jv[i]){
				tj = sum(tT[i,]<=TT[i,j],na.rm=TRUE)
				mu = xiv[tU[i,tj]]+as.vector(X[i,j,]%*%bev)
				si2 = si2+(Y[i,j]-mu)^2/sum(jv)
			}	
		}		
	}  
# Compute first log-likelihood component
	lk1 = sum(tnv*log(piv))
	for(u in 1:k) lk1 = lk1+sum(tN[u,-u]*log(Q[u,-u]))
	lk1 = lk1-sum(sv*qv)
# Compute second and third log-likelihood components  
	lk23 = lk_comp_joint_HM(par,si2,tT,tU,TT,W,X,Y,tv,dev,k,type)
	lk = lk1+lk23
	print("Final log-likelihood")
	print(lk)
# output
	out = list(bev=bev,nu=nu,xiv=xiv,phi=phi,psiv=psiv,piv=piv,Q=Q,si2=si2,tnv=tnv,tN=tN,sv=sv)
    
}
# ---------------------------------------------------------------------------------------------------
lk_comp_joint_HM <- function(par,si2,tT,tU,TT,W,X,Y,tv,dev,k,type){

# to compute the second and third components of the expect complete log-likelihood

# Preliminaries
	ncovW = ncol(W)
	ncovX = dim(X)[3]
	bev = par[1:ncovX]; eta = par[ncovX+1]; nu = exp(eta)
	xiv = par[(ncovX+2):(k+ncovX+1)]; phi = par[k+ncovX+2]; psiv = par[(k+ncovX+3):length(par)]
	n = nrow(TT)
	jv = rowSums(!is.na(Y))
	if(type=="norm") lsi2pi = log(si2)+log(2*pi)
# compute functions
	lk2 = lk3 = 0
	for(i in 1:n){
		tmp0 = as.vector(W[i,]%*%psiv)
		for(j in 1:jv[i]){
			tj = sum(tT[i,]<=TT[i,j],na.rm=TRUE)
			mu = xiv[tU[i,tj]]+as.vector(X[i,j,]%*%bev)
			if(type=="norm") lk2 = lk2-1/2*(lsi2pi+(Y[i,j]-mu)^2/si2)
			if(type=="bino"){
				mu = exp(mu)/(1+exp(mu))
				lk2 = lk2+(Y[i,j]*log(mu)+(1-Y[i,j])*log(1-mu))
			}
		}
		di = sum(tT[i,]<tv[i],na.rm=TRUE)
		if(tT[i,di]>0) for(j in 2:di) lk3 = lk3-(exp(xiv[tU[i,j-1]]*phi+tmp0)*(tT[i,j]^nu-tT[i,j-1]^nu))
		lk3 = lk3+dev[i]*(eta+(nu-1)*log(tv[i])+xiv[tU[i,di]]*phi+tmp0)-exp(xiv[tU[i,di]]*phi+tmp0)*
		                 (tv[i]^nu-tT[i,di]^nu)
	}
	lk23 = as.vector(lk2+lk3)

}
# ---------------------------------------------------------------------------------------------------
dlk_comp_joint_HM <- function(par,si2,tT,tU,TT,W,X,Y,tv,dev,k,type){

# Compute derivative of lk_comp_joint_HM  
  
# Preliminaries
  ncovW = ncol(W)
  ncovX = dim(X)[3]
  bev = par[1:ncovX]; eta = par[ncovX+1]; nu = exp(eta)
  xiv = par[(ncovX+2):(k+ncovX+1)]; phi = par[k+ncovX+2]; psiv = par[(k+ncovX+3):length(par)]
  n = nrow(TT)
  jv = rowSums(!is.na(Y))
  if(type=="norm") lsi2pi = log(si2)+log(2*pi)
# compute functions
	dlk2bev = rep(0,ncovX)
	dlk3eta = 0
	dlk2xiv = dlk3xiv = rep(0,k)
	dlk3phi = 0
	dlk3psiv = rep(0,ncovW)
	for(i in 1:n){
		tmp0 = as.vector(W[i,]%*%psiv)
		for(j in 1:jv[i]){
		  tj = sum(tT[i,]<=TT[i,j],na.rm=TRUE)
		  mu = xiv[tU[i,tj]]+as.vector(X[i,j,]%*%bev)
		  if(type=="norm"){
				dlk2bev = dlk2bev+X[i,j,]*(Y[i,j]-mu)/si2
				dlk2xiv[tU[i,tj]] = dlk2xiv[tU[i,tj]]+(Y[i,j]-mu)/si2
			}
			if(type=="bino"){
				mu = exp(mu)/(1+exp(mu))
				dlk2bev = dlk2bev+X[i,j,]*(Y[i,j]-mu)
				dlk2xiv[tU[i,tj]] = dlk2xiv[tU[i,tj]]+(Y[i,j]-mu)
			}
		}
		di = sum(tT[i,]<tv[i],na.rm=TRUE)
		if(tT[i,di]>0) for(j in 2:di){
			if(tT[i,j-1]>0) tmp = log(tT[i,j-1]) else tmp = 0
			tmp1 = exp(xiv[tU[i,j-1]]*phi+tmp0)
			dlk3eta = dlk3eta-tmp1*(log(tT[i,j])*tT[i,j]^nu-tmp*tT[i,j-1]^nu)*nu
			dlk3xiv[tU[i,j-1]] = dlk3xiv[tU[i,j-1]]-tmp1*(tT[i,j]^nu-tT[i,j-1]^nu)*phi
			dlk3phi = dlk3phi-tmp1*(tT[i,j]^nu-tT[i,j-1]^nu)*xiv[tU[i,j-1]]
			dlk3psiv = dlk3psiv-W[i,]*tmp1*(tT[i,j]^nu-tT[i,j-1]^nu)
		}
		if(tT[i,di]>0) tmp = log(tT[i,di]) else tmp = 0
		tmp1 = exp(xiv[tU[i,di]]*phi+tmp0)
		dlk3eta = dlk3eta+dev[i]*(1+nu*log(tv[i]))-tmp1*(log(tv[i])*tv[i]^nu-tmp*tT[i,di]^nu)*nu
		dlk3xiv[tU[i,di]] = dlk3xiv[tU[i,di]]+(dev[i]-tmp1*(tv[i]^nu-tT[i,di]^nu))*phi
		dlk3phi = dlk3phi+(dev[i]*xiv[tU[i,di]]-tmp1*(tv[i]^nu-tT[i,di]^nu)*xiv[tU[i,di]])
		dlk3psiv = dlk3psiv+W[i,]*(dev[i]-tmp1*(tv[i]^nu-tT[i,di]^nu))
	}
	dlk23 = c(dlk2bev,dlk3eta,dlk2xiv+dlk3xiv,dlk3phi,dlk3psiv)

}