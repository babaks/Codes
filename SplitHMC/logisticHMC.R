# This file contains all the files for running the standard HMC and the two Split HMC methods. This file is called by the R scripts for the examples discussed in the paper.  

# To run our codes, you need these two packages. Please make sure they are installed on your machine. 
library(mvtnorm)
library(numDeriv)


####### Standard HMC ##########


hmcStandard <- function(x, y, q.hat, n.iter, epsilon, L){

n.obs = dim(x)[1]	
n.param = dim(x)[2]	

post.samp <- matrix(NA, n.iter, n.param)

current.q <- q.hat
current.U <- U(x, y, current.q, sigma)
current.g = grad.U(x, y, current.q, sigma)

post.samp[1, ] <- current.q
for (i in 2:n.iter){
	samp <- getSamp(x, y, current.q, current.U, current.g, sigma, n.param, epsilon, L)
	current.q <- samp$q
	current.U <- samp$U
	current.g <- samp$g
	post.samp[i, ] <- current.q
	
	}
return(post.samp)
}

getSamp <- function(x, y, current.q, current.U, current.g, sigma, n.param, e, L){
	
q = current.q	
p = rnorm(n.param)
g = current.g

current.p = p
current.K = sum(current.p^2)/2
current.H = current.U+current.K

e = runif(1, 0.8*e, e)

for(leap in 1:L){
	p = p - e*g/2
    q = q + e*p
    g = grad.U(x, y, q, sigma)
    p = p - e*g/2
}
proposed.U = U(x, y, q, sigma)
proposed.K = sum(p^2)/2 
proposed.H = proposed.U+proposed.K

acceptProb = min(1, exp(current.H - proposed.H))
if(runif(1)<acceptProb){
    return (list(q = q, U = proposed.U, g = g))  # accept
}else{
    return (list(q = current.q, U = current.U, g = current.g))  # reject
}
}


U = function(x, y, beta, sigma){

eta = x%*%beta

logLike =  sum(y*eta - log(1+exp(eta)))

logPrior =  sum( (-0.5*(beta)^2 ) / (sigma^2) )

E = -(logLike + logPrior)

}

grad.U = function(x, y, beta, sigma){

eta = x%*%beta;

dLogLike = t(x)%*%(y - exp(eta - log(1 + exp(eta)) ))
dLogPrior =  -(beta) / (sigma^2)
 
g = -(dLogLike + dLogPrior);


}





####### Split HMC with Normal Approximation ##########


hmcSplitNorm <- function(x, y, q.hat, FI, R, n.iter, epsilon, L){

n.obs = dim(x)[1]	
n.param = dim(x)[2]	

post.samp <- matrix(NA, n.iter, n.param)

current.q <- q.hat
current.U <- U.norm(x, y, current.q, sigma)
current.g = grad.U1.norm(x, y, current.q, sigma, q.hat, FI) 

post.samp[1, ] <- current.q
for (i in 2:n.iter){
	samp <- getSampNorm(x, y, current.q, current.U, current.g, sigma, q.hat, FI, R, n.param, epsilon, L)
	current.q <- samp$q
	current.U <- samp$U
	current.g <- samp$g
	post.samp[i, ] <- current.q
	
	}
return(post.samp)
}



getSampNorm = function (x, y, current.q, current.U, current.g, sigma, q.hat, FI, R, n.param, e, L)
{	
	
    q = current.q
    p = rnorm(n.param)	
    g = current.g
  
    current.p = p
    current.K = sum(current.p^2)/2
    current.H = current.U+current.K
  
e = runif(1, 0.8*e, e)

  for (i in 1:L)
  {
  	# Outer splitting
    p = p - e * g / 2 # The first half step for momentum with potential U1
    # Inner splitting
    #X0 <- c(q - q.hat,p)

	#direct solution by diagonizing A
    X <-  R%*% c(q - q.hat,p)  # solution of dX/dt = A %*% X    
    q <- X[1:n.param]+q.hat
    p <- X[(n.param+1):(2*n.param)]  
    # End of inner spliting
    
    g = grad.U1.norm(x, y, q, sigma, q.hat, FI) 
    p = p - e * g/ 2 # The last half step for momentum with potential U1
  }
  proposed.U = U.norm(x, y, q, sigma)  
  proposed.K = sum(p^2)/2
  proposed.H = proposed.U + proposed.K

  acceptProb = min(1, exp(current.H - proposed.H))
  if (runif(1) < acceptProb)
  {
    return (list(q = q, U = proposed.U, g = g))  # accept
  }
  else
  {
    return (list(q = current.q, U = current.U, g = current.g))  # reject
  }
}



U.norm = function(x, y, beta, sigma){

eta = x%*%beta

logLike =  sum(y*eta - log(1+exp(eta)))

logPrior =  sum( (-0.5*(beta)^2 ) / (sigma^2) )

E = -(logLike + logPrior)

}


grad.U1.norm = function(x, y, beta, sigma, beta.hat, FI){

g0 = t(t(beta - beta.hat)%*%FI)

eta = x%*%beta;

dLogLike = t(x)%*%(y - exp(eta - log(1 + exp(eta)) ))
dLogPrior =  -(beta) / (sigma^2)

g = -(dLogLike + dLogPrior) - g0


}




####### Split HMC with Data Splitting ##########


hmcSplitData <- function(x, y, x0, y0, x1, y1, q.hat, n.iter, epsilon, L, M){

n.param = dim(x)[2]	

post.samp <- matrix(NA, n.iter, n.param)

current.q <- q.hat
post.samp[1, ] <- current.q

current.U = U.data(x, y, current.q, sigma) 

current.g = grad.U1.data(x1, y1, current.q, sigma)

for (i in 2:n.iter){
	samp <- getSampData(x, y, x0, y0, x1, y1, current.q, current.U, current.g, sigma, n.param, epsilon, L, M)
	current.q <- samp$q
	current.U <- samp$U
	current.g <- samp$g
	
	post.samp[i, ] <- current.q
	
	}
return(post.samp)
}



getSampData = function (x, y, x0, y0, x1, y1, current.q, current.U, current.g, sigma, n.param, e, L, M)
{
	
  q = current.q	
  g = current.g 
  p = rnorm(n.param)

  current.p = p
  
   
  current.K = sum(current.p^2)/2
  current.H = current.U+current.K
 
  e = runif(1, 0.8*e, e)
 
  newG1 = g
  for (i in 1:L)
  {
  	# Outer splitting

    p = p - e * newG1 / 2 # The first half step for momentum with potential U1

    # Inner splitting
    newG0 = grad.U0.data(x0, y0, q, sigma)
    for (m in 1:M){ # update with potential U0

    	p = p - e / (2 * M) * newG0
    	q = q + (e / M) * p
    	newG0 = grad.U0.data(x0, y0, q, sigma)
    	p = p - e / (2 * M) * newG0
    }

	newG1 = grad.U1.data(x1, y1, q, sigma)
    p = p - e * newG1 / 2 # The last half step for momentum with potential U1
  }

  # Evaluate potential and kinetic energies at start and end of trajectory

proposed.U = U.data(x, y, q, sigma) 

proposed.K = sum(p^2)/2
proposed.H = proposed.U + proposed.K

  acceptProb = min(1, exp(current.H - proposed.H))
  if (runif(1) < acceptProb)
  {
    return (list(q = q, U = proposed.U, g = newG1, Ind = 1))  # accept

  }
  else
  {
    return (list(q = current.q, U = current.U, g = current.g, Ind = 0))  # reject

  }
}



U.data = function(x, y, beta, sigma){

eta = x%*%beta

logLike =  sum(y*eta - log(1+exp(eta)))

logPrior =  sum( (-0.5*(beta)^2 ) / (sigma^2) )

E = -(logLike + logPrior)

}


grad.U0.data = function(x, y, beta, sigma){

eta = x%*%beta

dLogLike = t(x)%*%(y - exp(eta - log(1 + exp(eta)) ))
dLogPrior =  -(beta) / (sigma^2)
 
g = -(dLogLike + dLogPrior)

}



grad.U1.data = function(x, y, beta, sigma){


eta = x%*%beta

dLogLike = t(x)%*%(y - exp(eta - log(1 + exp(eta)) ))

g = -(dLogLike)

}





############### Measuring Performance ##################


# This function is used to get the autocorrelation time
ac.time <- function(x, q, m){

n = length(x)
s2 = var(x)
batch.mean = NULL
for(i in 1:q){
	batch.mean[i] <- mean(x[(((i-1)*m)+1):(i*m)])
	}	
return(m*var(batch.mean)/s2)

}

get.logLike = function(beta){

eta = x%*%beta

logLike =  sum(y*eta - log(1+exp(eta)))

}


# This function returns the performance measures
getPerformance <- function(samp, cpu.time, nIter){


logLike <- NULL
for(i in 1:dim(samp)[1]){

	logLike[i] <- get.logLike(samp[i, ])	
	
}

acp = length(unique(samp[, 1]))/length(samp[, 1])
CPU = round(cpu.time/nIter, 4)
batch.size = floor(nIter^(2/3))
n.batch = floor(nIter^(1/3))	
act = round(ac.time(logLike, n.batch, batch.size), 2)

efficiency = act*CPU
return(list(acp=acp, CPU=CPU, act = act, efficiency=efficiency))

}





######### MAP, Fisher Information, Splitting Data ###############



# This function is used to get the maximum a posteriori (MAP) estimate
get.map = function(beta){
	
eta = x%*%beta
logLike =  sum(y*eta - log(1+exp(eta)))
logPrior =  sum( (-0.5*(beta)^2 ) / (sigma^2) )
E = -(logLike + logPrior)

}

get.grad = function(beta){

eta = x%*%beta;

dLogLike = t(x)%*%(y - exp(eta - log(1 + exp(eta)) ))
dLogPrior =  -(beta) / (sigma^2)
 
g = -(dLogLike + dLogPrior);


}

get.hess = function(beta) {
		t(x)%*%diag(as.vector(exp(x%*%beta)/(1+exp(x%*%beta))^2))%*%x + diag(beta/(sigma^2))
	
	}


# Get MAP and Fisher Information
n.param = dim(x)[2]
opt.res <- optim(rnorm(n.param, 0, .2), get.map)
q.hat <- opt.res$par
FI <- get.hess(q.hat)
D = length(q.hat)
A <- matrix(0,2*D,2*D)
A[1:D,(D+1):(2*D)] <- diag(1,D)
A[(D+1):(2*D),1:D] <- -FI 
eig <- eigen(A)
lambda <- eig$values
Gamma <- eig$vectors
R = Re(Gamma %*% diag(exp(0.9*epsilon.splitNorm*lambda)) %*% solve(Gamma))


# Split the data
p =  (exp(x%*%q.hat)/(1+exp(x%*%q.hat)))
y.hat <- round(p)
p = cbind(p, 1-p)
ent <- -rowSums(p*log(p))
ent[is.na(ent)] <- 0
ind1 = which(ent > quantile(ent, frac))
U0.ind = ind1
U1.ind = setdiff(seq(1, dim(x)[1]), U0.ind)
x0 = x[U0.ind, ]
x1 = x[U1.ind, ]
y0 = y[U0.ind]
y1 = y[U1.ind]

