# system("R CMD SHLIB lk_rec.f --preclean")
# system("R CMD SHLIB lk_rec_surv.f --preclean")
# system("R CMD SHLIB lk_comp_joint.f --preclean")
# system("R CMD SHLIB lk_comp_joint_surv.f --preclean")
# system("R CMD SHLIB dlk_comp_joint.f --preclean")
# system("R CMD SHLIB dlk_comp_joint_surv.f --preclean")
library(expm)
library(survival)
library(numDeriv)
library(survival)

dyn.load("lk_rec.so")
dyn.load("lk_rec_surv.so")
dyn.load("lk_comp_joint.so")
dyn.load("lk_comp_joint_surv.so")
dyn.load("dlk_comp_joint.so")
dyn.load("dlk_comp_joint_surv.so")

est_jhmc_rec2fort = function(Y,TT,X,W,tv,censv,k,type="norm",btv=NULL,reltol=1e-6,inits=NULL,verbose=FALSE,
                             LC=FALSE,surv=TRUE,maxit=5000,stderr=FALSE) {

# fit joit model with continuos HM model
#
# INPUT:
#
# Y       = matrix of responses of dimensino n x Tmax (Tmax = maximum number of individual observations)
# TT      = corresponding matrix of observation times
# X       = array of covariates for longitudinal process of dimension n x TT x nX (nX = number of covariates X)
# W       = array of covariates (including itercept) for survival process of dimension n x nW (nW = number of covariates V)
# tv      = vector of surivival times of lenght n
# censv   = vector of indicators of censoring of lenght n
# k       = number of latent states
# type    = type of response variables (norm or bino)
# btv     = list of equally spaced points defining time windows for the likelihood approximation
# reltol  = relative tollerance for checking log-likelihood convergence
# inits   = list of initial values of the parameters
# verbose = to disply partial output
# LC      = to fit the latent class (instead of latent Markov) model
# surv    = if survival time is provided
#
# OUTPUT:
#
# lk      = maximum log-likelihood at convergence
# Q       = estimated infinitesimal transition matrix
# piv     = estimated vector of initial probabilities
# bev     = estiamted vector of regression coefficients of the longitudinal responses
# nu      = estiamted parater in the baseline hazard function
# xiv     = estimated vector of support points for the latent distribution
# phi     = estimated association coeffiecint for the survivial process
# psiv    = estimated vector of 
# si2     = estimated variance parameter
# pP      = time-specific posterior probabilities of each state
# Pi      = discrete time transition matrix
# it      = final number of iterations
# se      = standard errors (nu on the log scale)
# info    = information matrix (nu on the log scale)    

# Preliminaries
	t0 = proc.time()[3]
	npX = dim(X)[3]          # number of covariates on longitudinal process
	if(surv) npW = ncol(W)   # number of covariates on survival process
	n = dim(Y)[1]            # sample size
	jv = rowSums(!is.na(Y))  # number of longitudinal observations
	mjv = max(jv)            # maximum number of longitudinal observations
	if(!surv) tv = apply(TT,1,max,na.rm=TRUE)

# starting values
	Q = matrix(0,k,k)
	if(!is.null(inits)){
		if(LC) Q=matrix(0,k,k) else Q = inits$Q
		piv = inits$piv
		bev = inits$bev
		if(surv){
			nu = inits$nu
			phi = inits$phi
			psiv = inits$psiv
		}else{
			nu = phi = psiv = NULL
		}
		xiv = inits$xiv
		if(type=="bino") si2 = NULL else si2 = inits$si2
	}
	if(is.null(inits)){
		if(!LC){
			Q = matrix(1/(k-1),k,k)             
			diag(Q) = -1
		}
		if(type=="norm"){
			jnk = lm(Y[,1]~.,data=data.frame(X[,1,]))
			xiv = jnk$coef[1]+seq(-k,k,length=k)
			bev = jnk$coef[-1]
			res = jnk$residuals
			si2 = var(res)
		}
		if(type=="bino"){
			yv = as.vector(Y)
			XX = matrix(X,n*mjv,npX)
			jnk = glm(yv~XX,family=binomial)
			pred = predict(jnk)
			ppred0 = exp(-pred)/(1+exp(-pred))
			ppred1 = exp(pred)/(1+exp(pred))			
			res = yv[!is.na(yv)]*(-pred-log(1-ppred1)/ppred1)+(1-yv[!is.na(yv)])*(-pred+log(1-ppred0)/ppred0)
			si2 = NULL
		}
		out = kmeans(res,centers=k,iter.max=100)
		bev = jnk$coef[-1]
		ind = order(out$centers)
		xiv = jnk$coef[1]+out$centers[ind]
		piv = out$size[ind]/sum(out$size)
		if(surv){
			su = survreg(Surv(tv,ifelse(censv,0,1))~.,data=data.frame(W[,-1]))
			nu = 1/su$scale
			phi = 0
			jnk = su$coef
			jnk[-1] = jnk[-1]*nu
			jnk[1] = exp(jnk[1])
			psiv = -jnk
		}else{
			nu = phi = psiv = NULL			
		}
	}

# compute transition matrix
	if(is.null(btv)){
		a = min(c(0.1,tv))
		btv = seq(0,max(tv),a)
	}else{
		a = btv[2]-btv[1]
	}
	M = length(btv)
	if(LC){
		Pi = diag(k)
	}else{
		Pi = diag(k)+a*Q
		Tmp = a*Q
		it = 1
		while(any(abs(Tmp)>0) & it<1000){
			it = it+1
			Tmp = (a/it)*(Tmp%*%Q)
			Pi = Pi+Tmp
		}
	}

# expand matrix of responses and array of covariates
	YM = matrix(0,n,M)
	XM = array(0,c(n,M,npX))
	PM = matrix(0,n,M)
	for(i in 1:n) for(j in 1:jv[i]){
		m = sum(TT[i,j]>=btv)
		YM[i,m] = Y[i,j]
		XM[i,m,] = X[i,j,]
		PM[i,m] = 1 
	}
	
# compute log-likelihood
	if(is.null(si2)) si2f = as.double(0) else si2f = si2
	if(surv){
		out = .Fortran("lk_rec_surv",k=as.integer(k),piv=piv,Pi=Pi,npX=as.integer(npX),bev=bev,
		                nu=nu,xiv=xiv,phi=phi,npW=as.integer(npW),psiv=psiv,si2=si2f,
		                typ=type,n=as.integer(n),XM=XM,W=W,YM=YM,M=as.integer(M),PM=as.integer(PM),
		                btv=btv,dev=as.integer(1-censv),tv=as.double(tv),lk=0,Pp1=array(999,c(k,M,n)),Pp2=matrix(0,k,k))
	}else{
		out = .Fortran("lk_rec",k=as.integer(k),piv=piv,Pi=Pi,npX=as.integer(npX),bev=bev,
		                xiv=xiv,si2=si2f,typ=type,n=as.integer(n),XM=XM,YM=YM,M=as.integer(M),PM=as.integer(PM),
		                btv=btv,tv=as.double(tv),lk=0,Pp1=array(999,c(k,M,n)),Pp2=matrix(0,k,k),ver=as.integer(rep(0,n)))
	}
	lk = out$lk
	it = 0; lko = lk
	if(verbose) print(c(0,lk,phi))
	while(((lk-lko)/abs(lko)>reltol | it==0) & it<maxit){
		it = it+1; lko = lk
# M-step
# update piv and Pi
		piv = rowSums(out$Pp1[,1,])/n
		if(!LC) Pi = (1/rowSums(out$Pp2))*out$Pp2
# update other parameters
		if(surv) par = c(bev,log(nu),xiv,phi,psiv) else par = c(bev,xiv,phi,psiv)
		if(it==1){
			if(surv){
	 			par1 = optim(par,lk_comp_joint_HM2fort,method="Nelder-Mead",piv=piv,Pi=Pi,si2=si2,type=type,
	      		             dev=1-censv,tv=tv,XM=XM,W=W,YM=YM,PM=PM,btv=btv,Pp1=out$Pp1,
		     	             control=list(fnscale=-1))
		   	}else{
	 			par1 = optim(par,lk_comp_joint_HM2fort,method="Nelder-Mead",piv=piv,Pi=Pi,si2=si2,type=type,
	      		             tv=tv,XM=XM,YM=YM,PM=PM,btv=btv,Pp1=out$Pp1,surv=FALSE,
	      		             control=list(fnscale=-1))
		   	}
			par = par1$par
		}else{
			if(surv){
				par2 = optim(par,lk_comp_joint_HM2fort,dlk_comp_joint_HM2fort,method="BFGS",piv=piv,Pi=Pi,si2=si2,
	     		             type=type,dev=1-censv,tv=tv,XM=XM,W=W,YM=YM,PM=PM,btv=btv,Pp1=out$Pp1,
		     	             control=list(fnscale=-1,reltol=10^-12))
			}else{
				par2 = optim(par,lk_comp_joint_HM2fort,dlk_comp_joint_HM2fort,method="BFGS",piv=piv,Pi=Pi,si2=si2,
	     		             type=type,tv=tv,XM=XM,YM=YM,PM=PM,btv=btv,Pp1=out$Pp1,surv=FALSE,
		     	             control=list(fnscale=-1,reltol=10^-12))				
			}
			par = par2$par
		}
		bev = par[1:npX]
		if(surv){
			nu = exp(par[npX+1]); xiv = par[(npX+2):(k+npX+1)]
			phi = par[k+npX+2]; psiv = par[(k+npX+3):length(par)]
		}else{
			xiv = par[(npX+1):(k+npX)]
		}
		if(type=="norm"){
			num = 0
			for(i in 1:n) for(m in 1:M) if(PM[i,m]==1){
				muv = xiv+as.vector(XM[i,m,]%*%bev)
				num = num+(YM[i,m]-muv)^2%*%out$Pp1[,m,i]
			}
			si2 = as.vector(num/sum(jv))
			if(is.na(si2)){
				print(si2)
				browser()
			}
		}
# compute log-likelihood
		if(is.null(si2)) si2f = as.double(0) else si2f = si2
		if(surv){
			out = .Fortran("lk_rec_surv",k=as.integer(k),piv=piv,Pi=Pi,npX=as.integer(npX),bev=bev,
		                           nu=nu,xiv=xiv,phi=phi,npW=as.integer(npW),psiv=psiv,si2=si2f,
		                           typ=type,n=as.integer(n),XM=XM,W=W,YM=YM,M=as.integer(M),PM=as.integer(PM),
		                           btv=btv,dev=as.integer(1-censv),tv=as.double(tv),lk=0,Pp1=array(999,c(k,M,n)),Pp2=matrix(0,k,k))
		}else{
			out = .Fortran("lk_rec",k=as.integer(k),piv=piv,Pi=Pi,npX=as.integer(npX),bev=bev,
		                           xiv=xiv,si2=si2f,typ=type,n=as.integer(n),XM=XM,YM=YM,M=as.integer(M),PM=as.integer(PM),
		                           btv=btv,tv=as.double(tv),lk=0,Pp1=array(999,c(k,M,n)),Pp2=matrix(0,k,k))
		}
		lk = out$lk
		if(verbose & it%%100==0) print(c(it,lk,lk-lko,phi,proc.time()[3]-t0))
        }
	if(verbose) print(c(it,lk,lk-lko,phi,proc.time()[3]-t0))

# final output
if(!LC) {Q = logm(Pi %^% (1/a))}
    info=se=score=bread=meat=NULL
    if(stderr) {
        meat = Vdlk_comp_joint_HM2(par,piv,Pi,si2,type,XM,W,YM,jv,btv,out$Pp1,dev=as.integer(1-censv),tv=as.double(tv),surv)
        bread=-jacobian(function(x) {
           dlk_comp_joint_HM2fort(x,piv,Pi,si2,type,XM,W,YM,PM,jv,btv,out$Pp1,dev=as.integer(1-censv),tv=as.double(tv),surv)
},par)
       
        score=dlk_comp_joint_HM2fort(par,piv,Pi,si2,type,XM,W,YM,PM,jv,btv,out$Pp1,dev=as.integer(1-censv),tv=as.double(tv),surv)
        info = solve(bread)%*%meat%*%solve(bread)
	se=sqrt(diag(info))
	if(rcond(info)<1e-6) {
		se=pmax(sqrt(diag(solve(bread))),sqrt(diag(solve(meat))))}
    }
           
	out = list(lk=lk,Q=Q,piv=piv,bev=bev,nu=nu,xiv=xiv,phi=phi,psiv=psiv,si2=si2,pP=out$Pp1,Pi=Pi,it=it,tim=proc.time()[3]-t0,info=solve(info),se=se,score=score,bread=bread,meat=meat)
	return(out)

}
lk_comp_joint_HM2fort <- function(par,piv,Pi,si2,type,XM,W,YM,PM,btv,Pp1,dev,tv,surv=TRUE){

# to compute the second and third components of the expect complete log-likelihood

# preliminaries
	k = length(piv)
	npX = dim(XM)[3]
	if(surv) npW = ncol(W)
	bev = par[1:npX]
	if(surv){
		lnu = par[npX+1]; xiv = par[(npX+2):(k+npX+1)]; phi = par[k+npX+2]
		nu = exp(lnu); psiv = par[(k+npX+3):length(par)]
	}else{
		xiv = par[(npX+1):(k+npX)]
	}
	n = nrow(YM)
	M = length(btv)
# call fortran function
	if(is.null(si2)) si2f = as.double(0) else si2f = si2
	if(surv){
		out = .Fortran("lk_comp_joint_surv",bev=bev,lnu=lnu,xiv=xiv,phi=phi,psiv=psiv,k=as.integer(k),piv=piv,Pi=Pi,si2=si2f,
		                   typ=type,n=as.integer(n),M=as.integer(M),npX=as.integer(npX),XM=XM,npW=as.integer(npW),
	    	               W=W,YM=YM,PM=as.integer(PM),btv=btv,Pp1=Pp1,dev=as.integer(dev),tv=as.double(tv),lk=0)
	}else{
		out = .Fortran("lk_comp_joint",bev=bev,xiv=xiv,k=as.integer(k),piv=piv,Pi=Pi,si2=si2f,
		                   typ=type,n=as.integer(n),M=as.integer(M),npX=as.integer(npX),XM=XM,
	    	               YM=YM,PM=as.integer(PM),btv=btv,Pp1=Pp1,tv=as.double(tv),lk=0)
	}
	lk = out$lk
# return output	
	return(lk)

}
dlk_comp_joint_HM2fort <- function(par,piv,Pi,si2,type,XM,W,YM,PM,jv,btv,Pp1,dev,tv,surv=TRUE){

# to compute the derivative of lk_comp_joint_HM
# preliminaries
	k = length(piv)
	npX = dim(XM)[3]
	if(surv) npW = ncol(W)
	bev = par[1:npX]
	if(surv){
		lnu = par[npX+1]; xiv = par[(npX+2):(k+npX+1)]; phi = par[k+npX+2]
		nu = exp(lnu); psiv = par[(k+npX+3):length(par)]
	}else{
		xiv = par[(npX+1):(k+npX)]
	}
	n = nrow(YM)
	M = length(btv)
# call Fortran
# call fortran function
	if(is.null(si2)) si2f = as.double(0) else si2f = as.double(si2)
	if(surv){
		out = .Fortran("dlk_comp_joint_surv",bev=bev,lnu=lnu,xiv=xiv,phi=phi,psiv=psiv,k=as.integer(k),piv=piv,Pi=Pi,si2=si2f,
		                   typ=type,n=as.integer(n),M=as.integer(M),npX=as.integer(npX),XM=XM,npW=as.integer(npW),
	    	               W=W,YM=YM,PM=as.integer(PM),btv=btv,Pp1=Pp1,dev=as.integer(dev),tv=as.double(tv),
	    	               dlk=rep(0,npX+1+k+1+npW))
	}else{
		out = .Fortran("dlk_comp_joint",bev=bev,xiv=xiv,k=as.integer(k),piv=piv,Pi=Pi,si2=si2f,
		                   typ=type,n=as.integer(n),M=as.integer(M),npX=as.integer(npX),XM=XM,
	    	               YM=YM,PM=as.integer(PM),btv=btv,Pp1=Pp1,tv=as.double(tv),dlk=rep(0,npX+k))
	}
	dlk = out$dlk
# return output	
	return(dlk)


}



lk_comp_joint_HM2 <- function(par,piv,Pi,si2,type,XM,W,YM,btv,Pp1,dev,tv,surv=TRUE){

# to compute the second and third components of the expect complete log-likelihood

# preliminaries
	k = length(piv)
	npX = dim(XM)[3]
	bev = par[1:npX]
	if(surv){
		lnu = par[npX+1]; xiv = par[(npX+2):(k+npX+1)]; phi = par[k+npX+2]
		nu = exp(lnu); psiv = par[(k+npX+3):length(par)]
	}else{
		xiv = par[(npX+1):(k+npX)]
	}
	n = nrow(YM)
	M = length(btv)
	if(surv){
		btvnu = btv^nu
		dbtvnu = diff(btvnu)
	}
	if(type=="norm") lsi2pi = log(si2)+log(2*pi)
# compute functions
	l2 = 0
	if(surv){
		l3 = 0
		tmp00 = xiv*phi
	}
	if(surv) Wp = W%*%psiv
	for(i in 1:n){
		if(surv){
			tmp0 = tmp00+Wp[i]
			etmp0 = exp(tmp0)
		}
		mi = sum(btv<=tv[i])
		for(m in 1:mi){
			if(!is.na(YM[i,m])){
				muv = xiv+as.vector(XM[i,m,]%*%bev)
				if(type=="norm") l2 = l2-1/2*(lsi2pi+(YM[i,m]-muv)^2/si2)%*%Pp1[,m,i]
				if(type=="bino"){
					muv = exp(muv)/(1+exp(muv))
					l2 = l2+(YM[i,m]*log(muv)+(1-YM[i,m])*log(1-muv))%*%Pp1[,m,i]
				}
			}
		}
		if(surv){
			if(mi>1) for(m in 2:mi) l3 = l3-dbtvnu[m-1]*(etmp0%*%Pp1[,m-1,i])
			l3 = l3+(dev[i]*(lnu+(nu-1)*log(tv[i])+tmp0)-etmp0*(tv[i]^nu-btvnu[mi]))%*%Pp1[,mi,i]
		}
	}
	if(surv) lk = as.vector(l2+l3) else lk = l2
	return(lk)

}
dlk_comp_joint_HM2 <- function(par,piv,Pi,si2,type,XM,W,YM,jv,btv,Pp1,dev,tv,surv=TRUE){

# to compute the derivative of lk_comp_joint_HM
# preliminaries
	k = length(piv)
	npX = dim(XM)[3]
	if(surv) npW = ncol(W)
	bev = par[1:npX]
	if(surv){
		lnu = par[npX+1]; xiv = par[(npX+2):(k+npX+1)]; phi = par[k+npX+2]
		nu = exp(lnu); psiv = par[(k+npX+3):length(par)]
	}else{
		xiv = par[(npX+1):(k+npX)]
	}
	n = nrow(YM)
	M = length(btv)
	if(surv){
		btvnu = btv^nu
		dbtvnu = diff(btvnu)
	}
# compute functions
	dl2bev = rep(0,2)
	if(surv) dl3nu = 0
	dl2xiv = dl3xiv = rep(0,k)
	if(surv){
		dl3phi = 0
		dl3psiv = rep(0,npW)
	}
	for(i in 1:n){
		if(surv) tmp0 = exp(xiv*phi+as.vector(W[i,]%*%psiv))
		mi = sum(btv<=tv[i])
		for(m in which(!is.na(YM[i,1:mi]))){
			muv = xiv+as.vector(XM[i,m,]%*%bev)
			if(type=="norm"){
				dl2bev = dl2bev+XM[i,m,]*as.vector(((YM[i,m]-muv)/si2)%*%Pp1[,m,i])
				dl2xiv = dl2xiv+(YM[i,m]-muv)/si2*Pp1[,m,i]
			}
			if(type=="bino"){
				muv = exp(muv)/(1+exp(muv))
				tmp1 = (YM[i,m]-muv)*Pp1[,m,i]
				dl2bev = dl2bev+XM[i,m,]*sum(tmp1)
				dl2xiv = dl2xiv+tmp1
			}
		}
		if(surv){
			if(mi>1) for(m in 2:mi){
				if(btv[m-1]>0) tmp = log(btv[m-1]) else tmp = 0
				dl3nu = dl3nu-(tmp0*(log(btv[m])*btvnu[m]-tmp*btvnu[m-1]))%*%Pp1[,m-1,i]
				dl3xiv = dl3xiv-(tmp0*dbtvnu[m-1])*phi*Pp1[,m-1,i]
				dl3phi = dl3phi-(tmp0*dbtvnu[m-1]*xiv)%*%Pp1[,m-1,i]
				dl3psiv = dl3psiv-W[i,]*as.vector((tmp0*dbtvnu[m-1])%*%Pp1[,m-1,i])
			}
			if(btv[mi]>0) tmp = log(btv[mi]) else tmp = 0          
			dl3nu = dl3nu+(dev[i]*(1/nu+log(tv[i]))-tmp0*(log(tv[i])*tv[i]^nu-tmp*btvnu[mi]))%*%Pp1[,mi,i]
			dl3xiv = dl3xiv+(dev[i]*phi-tmp0*(tv[i]^nu-btvnu[mi])*phi)*Pp1[,mi,i]
			dl3phi = dl3phi+(dev[i]*xiv-tmp0*(tv[i]^nu-btvnu[mi])*xiv)%*%Pp1[,mi,i]
			dl3psiv = dl3psiv+W[i,]*as.vector((dev[i]-tmp0*(tv[i]^nu-btvnu[mi]))%*%Pp1[,mi,i])
		}
	}
# output
	if(surv) dlk = c(dl2bev,dl3nu*nu,dl2xiv+dl3xiv,dl3phi,dl3psiv) else dlk = c(dl2bev,dl2xiv)
	return(dlk)

}

Vdlk_comp_joint_HM2 <- function(par,piv,Pi,si2,type,XM,W,YM,jv,btv,Pp1,dev,tv,surv=TRUE){

# to compute the derivative of lk_comp_joint_HM
# preliminaries
	k = length(piv)
	npX = dim(XM)[3]
	if(surv) npW = ncol(W)
	bev = par[1:npX]
	if(surv){
		lnu = par[npX+1]; xiv = par[(npX+2):(k+npX+1)]; phi = par[k+npX+2]
		nu = exp(lnu); psiv = par[(k+npX+3):length(par)]
	}else{
		xiv = par[(npX+1):(k+npX)]
	}
	n = nrow(YM)
	M = length(btv)
	if(surv){
		btvnu = btv^nu
		dbtvnu = diff(btvnu)
	}
# compute functions
	
    for(i in 1:n){
        dl2bev = rep(0,npX)
	if(surv) dl3nu = 0
	dl2xiv = dl3xiv = rep(0,k)
	if(surv){
		dl3phi = 0
		dl3psiv = rep(0,npW)
	}
		if(surv) tmp0 = exp(xiv*phi+as.vector(W[i,]%*%psiv))
		mi = sum(btv<=tv[i])
		for(m in which(!is.na(YM[i,1:mi]))){
			muv = xiv+as.vector(XM[i,m,]%*%bev)
			if(type=="norm"){
				dl2bev = dl2bev+XM[i,m,]*as.vector(((YM[i,m]-muv)/si2)%*%Pp1[,m,i])
				dl2xiv = dl2xiv+(YM[i,m]-muv)/si2*Pp1[,m,i]
			}
			if(type=="bino"){
				muv = exp(muv)/(1+exp(muv))
				tmp1 = (YM[i,m]-muv)*Pp1[,m,i]
				dl2bev = dl2bev+XM[i,m,]*sum(tmp1)
				dl2xiv = dl2xiv+tmp1
			}
		}
		if(surv){
			if(mi>1) for(m in 2:mi){
				if(btv[m-1]>0) tmp = log(btv[m-1]) else tmp = 0
				dl3nu = dl3nu-(tmp0*(log(btv[m])*btvnu[m]-tmp*btvnu[m-1]))%*%Pp1[,m-1,i]
				dl3xiv = dl3xiv-(tmp0*dbtvnu[m-1])*phi*Pp1[,m-1,i]
				dl3phi = dl3phi-(tmp0*dbtvnu[m-1]*xiv)%*%Pp1[,m-1,i]
				dl3psiv = dl3psiv-W[i,]*as.vector((tmp0*dbtvnu[m-1])%*%Pp1[,m-1,i])
			}
			if(btv[mi]>0) tmp = log(btv[mi]) else tmp = 0          
			dl3nu = dl3nu+(dev[i]*(1/nu+log(tv[i]))-tmp0*(log(tv[i])*tv[i]^nu-tmp*btvnu[mi]))%*%Pp1[,mi,i]
			dl3xiv = dl3xiv+(dev[i]*phi-tmp0*(tv[i]^nu-btvnu[mi])*phi)*Pp1[,mi,i]
			dl3phi = dl3phi+(dev[i]*xiv-tmp0*(tv[i]^nu-btvnu[mi])*xiv)%*%Pp1[,mi,i]
			dl3psiv = dl3psiv+W[i,]*as.vector((dev[i]-tmp0*(tv[i]^nu-btvnu[mi]))%*%Pp1[,mi,i])
		}

        if(surv) dlk = c(dl2bev,dl3nu*nu,dl2xiv+dl3xiv,dl3phi,dl3psiv) else dlk = c(dl2bev,dl2xiv)
        if(i==1) res = dlk%*%t(dlk)
        if(i>1) res = res + dlk%*%t(dlk)
        
	}
# output
	
	return(res)

}


lk_rec_joint_HM2 <- function(piv,Pi,bev,nu,xiv,phi,psiv,si2,type,XM,W,YM,btv,dev,tv,surv=TRUE){
	
# preliminaries
	n = nrow(YM)
	k = length(xiv)
	M = length(btv)
	if(type=="norm") si = sqrt(si2)
	if(surv){
		btvnu = btv^nu
		dbtvnu = diff(btvnu)
	}
# recursion for each individual
	lk = 0
	Pp1 = array(NA,c(k,M,n))
	Pp2 = matrix(0,k,k)
	for(i in 1:n){
		mi = sum(btv<=tv[i])
		if(surv) tmp0 = exp(xiv*phi+as.vector(W[i,]%*%psiv))
# FORWARD RECURSION	
		Frec = matrix(0,k,mi)
# first time
		muv = xiv+as.vector(XM[i,1,]%*%bev)
		if(type=="norm") frec = piv*dnorm(YM[i,1],muv,si)
		if(type=="bino"){
			muv = exp(muv)/(1+exp(muv))  
			frec = piv*(YM[i,1]*muv+(1-YM[i,1])*(1-muv))
		}
		Frec[,1] = frec
# following time occasions
		if(mi>1) for(m in 2:mi){
			if(surv){
				tmp = exp(-tmp0*dbtvnu[m-1])
				frec = as.vector(t(Pi)%*%(frec*tmp))
			}else{
				frec = as.vector(t(Pi)%*%frec)				
			}
			if(!is.na(YM[i,m])){
				muv = xiv+as.vector(XM[i,m,]%*%bev)
				if(type=="norm") frec = frec*dnorm(YM[i,m],muv,si)
				if(type=="bino"){
					muv = exp(muv)/(1+exp(muv))  
					frec = frec*(YM[i,m]*muv+(1-YM[i,m])*(1-muv))
				}
			}
			Frec[,m] = frec	
		}
		if(surv){
			tmp = exp(-tmp0*(tv[i]^nu-btvnu[mi]))
			if(dev[i]==1) tmp = tmp*nu*tv[i]^(nu-1)*tmp0	
			frec = as.vector(frec*tmp)
		}
		fi = sum(frec)
# accumulate log-likelihood
		lk = lk+log(fi)

# BACKWARD RECURSION
		Grec = matrix(0,k,mi)
# last time occasion
		if(surv){
			grec = exp(-tmp0*(tv[i]^nu-btvnu[mi]))
			if(dev[i]==1) grec = grec*nu*tv[i]^(nu-1)*tmp0
		}else{
			grec = rep(1,k)
		}
		Grec[,mi] = grec
# previous time occasions
		if(mi>1) for(m in (mi-1):1){
			if(surv) tmp = exp(-tmp0*dbtvnu[m])
			if(!is.na(YM[i,m+1])){
				muv = xiv+as.vector(XM[i,m+1,]%*%bev)
				if(type=="norm") grec = grec*dnorm(YM[i,m+1],muv,si)
				if(type=="bino"){
					muv = exp(muv)/(1+exp(muv))  
					grec = grec*(YM[i,m+1]*muv+(1-YM[i,m+1])*(1-muv))
				}
			}
			if(surv) grec = tmp*as.vector(Pi%*%grec) else grec = as.vector(Pi%*%grec)
			Grec[,m] = grec	
		}
# Posterior probabilities for the first time occasion
		Pp1[,1:mi,i] = Frec*Grec/fi
# Joint Posterior probabilities for two occasions
		if(mi>1) for(m in 1:(mi-1)){
			if(surv) tmp = exp(-tmp0*dbtvnu[m])
			if(!is.na(YM[i,m+1])){
				muv = xiv+as.vector(XM[i,m+1,]%*%bev)
				if(type=="norm") tmp1 = dnorm(YM[i,m+1],muv,si)
				if(type=="bino"){
					muv = exp(muv)/(1+exp(muv))  
					tmp1 = YM[i,m+1]*muv+(1-YM[i,m+1])*(1-muv)
				}
			}else{
				tmp1 = 1
			}
			if(surv){
				Tmp = ((Frec[,m]*tmp)%o%(Grec[,m+1]*tmp1))*Pi/fi
			}else{
				Tmp = (Frec[,m]%o%(Grec[,m+1]*tmp1))*Pi/fi
			}
			Pp2 = Pp2+Tmp
		}
	}
# output
	out = list(lk=lk,Pp1=Pp1,Pp2=Pp2)

}



