library(mvtnorm)

n = 50
set.seed(2)
x = runif(n, -2.5, 2.5)
#y = (1 + 2*x - 3*sin(x)) + rnorm(n, 0, .5)
y = .5*(x^3) - 3*x + 1 + rnorm(n, 0, .5)
plot(x, y, xlim=c(-3, 3))

nIter = 2000
burnIn = 500

gpReg = function(x, y, nIter){


diffMatAll = matrix(x, nrow=n, ncol=n) - matrix(x, nrow=n, ncol=n, byrow=TRUE)
diffMatAll2 = diffMatAll^2

post.lambda = rep(NA, nIter)
post.eta = rep(NA, nIter)
post.rho = rep(NA, nIter)
post.sigma = rep(NA, nIter)


lambda = 1
eta = 1
rho = 1
sigma = 1

for(i in 1:nIter){

	lambda <- get.lambda(x, y, diffMatAll2, lambda, eta, rho, sigma)

	eta <- get.eta(x, y, diffMatAll2, lambda, eta, rho, sigma)
	
	rho <- get.rho(x, y, diffMatAll2, lambda, eta, rho, sigma)

	sigma <- get.sigma(x, y, diffMatAll2, lambda, eta, rho, sigma)

	post.lambda[i] <- lambda	
	post.eta[i] <- eta
	post.rho[i] <- rho
	post.sigma[i] <- sigma

	}


return(list(lambda = post.lambda, eta = post.eta, rho = post.rho, sigma = post.sigma))

}



get.lambda <- function(x, y, diffMatAll2, lambda, eta, rho, sigma){

	w = 4
	m = 10
	
	z = getPost(x, y, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = lambda - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J
	
	L = max(0, L)
	while (J>0 && L>0 && z < getPost(x, y, diffMatAll2, L, eta, rho, sigma)) {
		L = L - w	
		L = max(0, L)		
		J = J - 1
	}

	while (K>0 && z < getPost(x, y, diffMatAll2, R, eta, rho, sigma)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)

	while (z > getPost(x, y, diffMatAll2, newParam, eta, rho, sigma)) {
		if (newParam < lambda) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}
	

	return(newParam)
}




get.eta <- function(x, y, diffMatAll2, lambda, eta, rho, sigma){

	w = 4
	m = 10
	
	z = getPost(x, y, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = eta - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J
	
	L = max(0, L)
	while (J>0 && L>0 && z < getPost(x, y, diffMatAll2, lambda, L, rho, sigma)) {
		L = L - w	
		L = max(0, L)		
		J = J - 1
	}

	while (K>0 && z < getPost(x, y, diffMatAll2, lambda, R, rho, sigma)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)
	
	while (z > getPost(x, y, diffMatAll2, lambda, newParam, rho, sigma)) {
		if (newParam < eta) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}

	return(newParam)
}







get.rho <- function(x, y, diffMatAll2, lambda, eta, rho, sigma){

	w = 4
	m = 10
	
	z = getPost(x, y, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = rho - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J
	
	L = max(0, L)
	while (J>0 && L>0 && z < getPost(x, y, diffMatAll2, lambda, eta, L, sigma)) {
		L = L - w	
		L = max(0, L)		
		J = J - 1
	}

	while (K>0 && z < getPost(x, y, diffMatAll2, lambda, eta, R, sigma)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)
	
	while (z > getPost(x, y, diffMatAll2, lambda, eta, newParam, sigma)) {
		if (newParam < rho) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}

	return(newParam)
}





get.sigma <- function(x, y, diffMatAll2, lambda, eta, rho, sigma){

	w = 4
	m = 10
	
	z = getPost(x, y, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = sigma - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J
	
	L = max(0, L)
	while (J>0 && L>0 && z < getPost(x, y, diffMatAll2, lambda, eta, rho, L)) {
		L = L - w	
		L = max(0, L)		
		J = J - 1
	}

	while (K>0 && z < getPost(x, y, diffMatAll2, lambda, eta, rho, R)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)
	
	while (z > getPost(x, y, diffMatAll2, lambda, eta, rho, newParam)) {
		if (newParam < sigma) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}

	return(newParam)
}




getPost = function(x, y, diffMatAll2, lambda, eta, rho, sigma){

C = lambda + eta*exp(-rho*(diffMatAll2)) + sigma*diag(1, nrow=n, ncol=n);

L = chol(C);
invL = solve(L);
invC = invL%*%t(invL)
invL.y = invL%*%y

logLike =  - 0.5*log(det(C)) - 0.5 * t(invL.y)%*%invL.y

logPrior =  dgamma(lambda, 1, 1, log=TRUE) + dgamma(eta, 1, 1, log=TRUE) + dgamma(rho, 1, 1, log=TRUE) + dgamma(sigma, 1, 1, log=TRUE)


logPost = (logLike + logPrior)

}




samp = gpReg(x, y, nIter)





pred = function(x.tr, y.tr, x.te, samp, burnIn, nIter){

lambda = samp$lambda[burnIn:nIter]
eta = samp$eta[burnIn:nIter]
rho = samp$rho[burnIn:nIter]

S = length(lambda)

x = c(x.tr, x.te);
nTrain = length(x.tr)
nTest = length(x.te)
n = length(x);


diffMatAll = matrix(x, nrow=n, ncol=n) - matrix(x, nrow=n, ncol=n, byrow=TRUE)

y.hat = NULL
v.hat = NULL
for(i in 1:S){
	
C = lambda[i] + eta[i]*exp(-rho[i]*(diffMatAll^2))+ sigma*diag(1, nrow=n, ncol=n);

Ctrn = C[1:nTrain, 1:nTrain];
invCtrn = solve(Ctrn)

K = C[1:nTrain, (nTrain+1):n];
v = C[(nTrain+1):n, (nTrain+1):n];

s = sqrt(diag(v))

# E(y.te | y.tr)
for(rep in 1:5){
y.hat = cbind(y.hat, rnorm(nTest, t(K)%*%invCtrn%*%y.tr, s))
}

}

return(y.hat)

}

plot(x, y, xlim=c(-3, 3), ylim=c(-3, 5))

x.te <- seq(-3, 3, .1)
res = pred(x, y, x.te, samp, burnIn, nIter)
s = apply(res, 1, sd)
y.hat = rowMeans(res)

lines(x.te, y.hat)

lines(x.te, y.hat + 2*s, lty=2)
lines(x.te, y.hat - 2*s, lty=2)

