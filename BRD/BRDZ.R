

BRD = function(beta, nIter=2000, burnIn=500) {
	
p = length(beta)

#Initial values for MCMC and fixed parameters
lambda = 0.5	# lambda is P(H_1)

# lambda ~ Beta(aBeta, bBeta)			
aBeta = 1
bBeta = 1

nu = -3
sc = 2

aGamma = -3
bGamma = 2
alpha.dp = 1


nComp = min(3, p)
J = sample(seq(1, nComp), p, replace=TRUE)
nj = as.vector(table(J))
		
tau2.star = exp(rnorm(nComp, nu, sc))

mu0 = 0 
tau2 = tau2.star[J]

# These are for storing posterior samples 
#post.tau2 = matrix(NA, nrow=nIter, ncol = p)
post.tau2.0 = rep(NA, nIter)
#post.J = matrix(0, nrow=nIter, ncol = p)
post.lambda = rep(NA, nIter)

post.rank = matrix(NA, nrow=nIter, ncol = p)
	
samp = list()


#MCMC samples
for (iter in 1:nIter){
 	#print(iter)
	mcmc.res = dpMCMC(beta, tau2.star, J, nj, alpha.dp, mu0, nu, sc) 
		
	J = mcmc.res$J
	tau2.star = mcmc.res$tau2.star
	nj = mcmc.res$nj
	nComp = length(nj)
				
	tau2.star = remix(beta, J, tau2.star, mu0, nu, sc)

 
 	post.rank[iter, ] = rank(tau2.star)[J]
	tau2 = tau2.star[J]

	alpha.dp = pickAlpha(alpha.dp, nComp, p, aGamma, bGamma)
 	
 	mat.nj = matrix(nj, nComp, p)
 	mat.tau2 = matrix(tau2.star, nComp, p)
 	mat.beta = matrix(beta, nComp, p, byrow=TRUE)
 	
		
	#post.tau2[iter, ] = tau2	
	post.lambda[iter] = mean(tau2==min(tau2.star))
	#post.J[iter, tau2==min(tau2.star)] = 1
	post.tau2.0[iter] = min(tau2.star)
}	

r.tab = apply(post.rank, 2, table)		
r = as.numeric(names(unlist(lapply(r.tab, which.max))))

	
return(list(post.lambda = post.lambda[burnIn:nIter], post.rank = post.rank, r = r, post.tau2.0 = post.tau2.0 ))
	
}



dpMCMC = function(beta.scaled, tau2.star, J, nj, alpha.dp, mu0, nu, sc){
M=5 # Number of auxillary components
p = length(beta.scaled)
p.H0 = NULL
p.H1 = NULL
for (i in 1:p){
	beta = beta.scaled[i]
	phi = NULL
    curInd = J[i]     	
    nj[curInd] = nj[curInd] - 1
    if (nj[curInd] == 0) {
    	phi = tau2.star[curInd]
        nj = nj[-curInd]
		tau2.star = tau2.star[-curInd]                        
        J[J>curInd] = J[J>curInd] - 1
        kBar = length(nj)
        tau2.star[kBar+1] = phi
        for (m in 1:(M-1)){
            tau2.star[kBar+1+m] = exp(rnorm(1, nu, sc))
        }
    }else{
        kBar = length(nj)
        for (m in 1:M) {
           	tau2.star[kBar+m] = exp(rnorm(1, nu, sc))
    	}
    }			
    q1 = rep(0, kBar) 
    for (k in 1:kBar){
        q1[k] = getLogLike(beta, mu0, tau2.star[k])      
    }
    q1 = q1 + log(nj) - log(p-1+alpha.dp)
    q2 = rep(0, M)
    for (k in 1:M){
        q2[k] = getLogLike(beta, mu0, tau2.star[kBar+k])
    }
    q2 = q2 + (log(alpha.dp) - log(M)) - log(p-1+alpha.dp)
    q = c(q1, q2)
    qMax = max(q)
    qRel = q - qMax
    q = exp(qRel)   
    q = q/sum(q)    

    picked = which(rmultinom(1, 1, q)==1)
    
    if (picked <= kBar){
        J[i] = picked
        nj[picked] = nj[picked]+1
        tau2.star = tau2.star[-((kBar+1):(kBar+M))]
        q = q[-((kBar+1):(kBar+M))]
		}else{	
        J[i] = kBar+1
        nj = c(nj, 1)   
       	phi = tau2.star[picked]
        tau2.star = tau2.star[-((kBar+1):(kBar+M))]
        q = q[-((kBar+1):(kBar+M))]
        tau2.star[kBar+1] = phi      
        }


}    
newNj = nj
newJ = J
return(list(tau2.star = tau2.star, J = J, nj = nj, p.H0 = p.H0, p.H1 = p.H1))
}# end dpMCMC	
	
	
getLogLike = function(beta, mu, tau2){

logLike = sum(dnorm(beta, mu, sqrt(tau2), log=TRUE))

return(logLike )

}


remix = function(beta.scaled, J, tau2.star, mu0, nu, sc){
tau2.mat = tau2.star[J]


for(i in 1:length(tau2.star)){
	beta = beta.scaled[J == i]
	tau2.star[i] = getTau2(beta, mu0, tau2.star[i], nu, sc)
	}
	
return(tau2.star)

}	
	
	
getTau2 = function(beta, mu, curTau2, nu, sc){
	w = 1
	m = 40
	ltau2 = log(curTau2)	
	z = getPostTau2(beta, mu, ltau2, nu, sc) - rexp(1)

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = ltau2 - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J

	while (J>0  && z < getPostTau2(beta, mu, L, nu, sc)) {
		L = L - w
		J = J - 1
	}

	while (K>0 && z < getPostTau2(beta, mu, R, nu, sc)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)
	
	while (z > getPostTau2(beta, mu, newParam, nu, sc)) {
		if (newParam < ltau2) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}

	return(exp(newParam))
	}


getPostTau2 = function(beta, mu, ltau2, nu, sc) {
	tau2 = exp(ltau2)
	logPrior = dnorm(ltau2, nu, sc, log=TRUE)
	logLike = sum(dnorm(beta, mu, sqrt(tau2), log=TRUE))
	return(logPrior+logLike)
}


pickAlpha = function(curAlpha, iStar, p, a0, b0){
	w = 2
	m = 20
	
	currentParam = log(curAlpha)
	z = getPostAlpha(currentParam, iStar, p, a0, b0) - rexp(1)
	

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = currentParam - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J

	while (J>0 && z < getPostAlpha(L, iStar, p, a0, b0)) {
		L = L - w
		J = J - 1
	}

	while (K>0 && z < getPostAlpha(R, iStar, p, a0, b0)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)
	
	while (z > getPostAlpha(newParam, iStar, p, a0, b0)) {
		if (newParam < currentParam) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}

	return(exp(newParam))
	}


getPostAlpha = function(lalpha, iStar, p, a0, b0) {

alpha = exp(lalpha)

logLike = iStar*log(alpha) +  lgamma(alpha) - lgamma(alpha+p)  

logPrior = log(dnorm(lalpha, a0, b0))

logPost = logLike + logPrior

return(logPost)
}


# to run the code, enter the following commands, where beta is the list of z-scores, nIter is the number of iterations, and burnIn is the number of burn-ins. 

# simRes is list of output objects. The main output is "post.rank" which has the posterior samples of gene ranks. Genes with higher degree of relevance to the outcome variable have higher posterior ranks. 

# simRes = BRD(z, nIter, burnIn)

