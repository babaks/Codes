


# This is the main part of the code for the NARD model.

BRD = function(y, c, nIter=2000, burnIn=200) {
	
# c is the class indicator, i.e., treated/non treated 
# y are gene expression values 
 
n = dim(y)[1]
p = dim(y)[2]

n1 = sum(c==1)
n2 = sum(c==2)

# Adding a column of 1's for the intercept
x = cbind(rep(1, n), c-1); 

# d is the number of regression parameters including the intercept.
d = dim(x)[2]

X = t(x)%*%x;
L.x = chol(X);
invX = chol2inv(L.x);
L.x.inv = chol(invX);
 
# This give the mean of the posterior distribution.
BETA.hat = invX%*%t(x)%*%y;

x.a = matrix(x[, 1:(d-1)], n, d-1, byrow=TRUE)
x.b = matrix(x[, d], n, 1)

X.a = t(x.a)%*%x.a;
L.a = chol(X.a);
invX.a = chol2inv(L.a);
L.a.inv = chol(invX.a);


X.b = t(x.b)%*%x.b;
L.b = chol(X.b);
invX.b = chol2inv(L.b);
L.b.inv = chol(invX.b);


# The following two matrices hold the posterior samples for sigma2 and beta.
alpha = matrix(BETA.hat[1:(d-1), ], nrow=d-1, ncol=p, byrow=TRUE)
beta = matrix(BETA.hat[d, ], 1, p)
sigma2 = colSums((y - x%*%BETA.hat)^2)/(n-1)

#Initial values for MCMC and fixed parameters
lambda = 0.5	# lambda is P(H_1)

# lambda ~ Beta(aBeta, bBeta)			
aBeta = 1
bBeta = 9

nu = -2
sc = 2

aGamma = -3
bGamma = 2
alpha.dp = 1

tau2 = rep(NA, p)

nComp = min(3, p)
J = sample(seq(1, nComp), p, replace=TRUE)
nj = as.vector(table(J))
		
tau2.star = exp(rnorm(nComp, nu, sc))
	
mu0 = 0 
tau2 = tau2.star[J]

# These are for storing posterior samples 
post.tau2 = matrix(NA, nrow=nIter, ncol = p)
post.tau2.0 = rep(NA, nIter)
post.J = matrix(0, nrow=nIter, ncol = p)
post.lambda = rep(NA, nIter)

p.val =  matrix(NA, nrow=nIter, ncol = p)
post.p.H0 = matrix(NA, nrow=nIter, ncol = p)
post.p.H1 =  matrix(NA, nrow=nIter, ncol = p)
post.rank = matrix(NA, nrow=nIter, ncol = p)
post.beta = matrix(NA, nIter, p)
	
samp = list()


#MCMC samples
for (iter in 1:nIter){
	#print(iter)
	y.a = y - x.b%*%beta

	alpha.hat = invX.a%*%t(x.a)%*%y.a

	sigma2.mat = matrix(sigma2, d-1, p, byrow=TRUE)
        
    u = matrix(rnorm((d-1)*p), p, d-1);
    alpha = (t(u%*%L.a.inv) + alpha.hat/sqrt(sigma2.mat))*sqrt(sigma2.mat)
        
	y.b = y - x.a%*%alpha
	y.b = y.b[c==2, ]
	
	V = 1/(1/(tau2*sigma2) + n2/sigma2)
	mu_n = V*(mu0/(tau2*sigma2) + colSums(y.b)/sigma2) 
	sigma_n = sqrt(V)
	beta = rnorm(p, mean = mu_n, sd = sigma_n)
	
	BETA = rbind(alpha, beta)

    eps = y - x%*%BETA
    eps.bar = colMeans(eps)
    
    nu_n = n-1
    sigma02_n = colSums( (eps - eps.bar)^2 ) /(n-1)
    z = rchisq(p, nu_n);
    sigma2 = nu_n*sigma02_n/z;

	beta.scaled = beta/sqrt(sigma2)
		    	 	
	mcmc.res = dpMCMC(beta.scaled, tau2.star, J, nj, alpha.dp, mu0, nu, sc) 
		
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
 	
	q = log(mat.nj) - log(p) + dnorm(mat.beta/sqrt(mat.tau2), log=TRUE)
	max.q = apply(q, 2, max)
	A = max.q + log(colSums(exp(q - matrix(max.q, nComp, p, byrow=TRUE))))
	probs = exp(q - matrix(A, nComp, p, byrow=TRUE))
	
 	post.p.H0[iter, ] = probs[which(tau2.star==min(tau2.star)), ]
 	
 	if(length(nj)>1){
 	post.p.H1[iter, ] = probs[which.max(tau2.star), ]
 	}else{
 	post.p.H1[iter, ] = 0
 	}
		
 	post.beta[iter, ] = beta
	post.tau2[iter, ] = tau2	
	post.lambda[iter] = mean(tau2==min(tau2.star))
	post.J[iter, tau2==min(tau2.star)] = 1
	post.tau2.0[iter] = min(tau2.star)			
	tau2.0 <- min(tau2.star)
	p.val[iter, ] = pchisq((beta.scaled/sqrt(tau2.0))^2, 1, lower.tail=FALSE)

}	

r.tab = apply(post.rank, 2, table)		
r = as.numeric(names(unlist(lapply(r.tab, which.max))))
	
return(list(post.lambda = mean(post.lambda[burnIn:nIter]), post.tau2.0 = post.tau2.0, post.tau2 = colMeans(post.tau2[burnIn:nIter, ]), post.beta = colMeans(post.beta[burnIn:nIter, ]), post.rank = post.rank, r = r ))
	
}


# This is the MCMC part for the Dirichlet process model

dpMCMC = function(beta.scaled, tau2.star, J, nj, alpha.dp, mu0, nu, sc){
M=3 # Number of auxillary components
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
