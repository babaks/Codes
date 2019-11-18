library(coda)


# Simulating the data: I set, lambda = 0.7, sigma2 = 0.5, eta2 = 8
x = c(rnorm(700, 0, sqrt(0.5)), rnorm(300, 0, sqrt(0.5 + 8)))
x = sample(x)
write.table(round(x, 2), file='tests.txt', row.names=FALSE, col.names=FALSE)



# These are the parameters of the prior for sigma2 and eta2
alpha1 = 1; beta1 = 0.5;

# These are the parameters of the prior for lambda
alpha2 = 2; beta2 = 8;


sigma = 3;
# This is the observed dataset
yObs = scan("tests.txt", what = numeric());
n = length(yObs);


assignment2Slice = function(nIterations) {

# These are the parameters of the prior for mu


# These are our initial points
sigma2 = numeric();
eta2 = numeric();
lambda = numeric();

sigma2[1] = 1
eta2[1] = 1;
lambda[1] = 0.5;

for (i in 2:nIterations) {

	#  As described by Neal (2003), "In practice, it is often safer to compute g(x) = log(f (x)) rather than f(x)  itself, in order to avoid possible problems with floating-point underflow. One can then use the auxiliary variable z = log(y) = g(x) -e, where e is exponentially distributed with mean one, and define the slice by S = {x : z < g(x)}."   


# Sampling sigma2
	
	w = 1; m = 20;
	 
	z = getPostSigma2(sigma2[i-1], eta2[i-1], lambda[i-1], yObs) - rexp(1);
	currentParam = sigma2[i-1];

	# Stepping out to obtain the [L, R] range
	u = runif(1);
	L = currentParam - w*u;
	R = L + w;
	v = runif(1);
	J = floor(m*v);
	K = (m-1) - J;

	L = max(0, L);
	while (J>0 && z < getPostSigma2(L, eta2[i-1], lambda[i-1], yObs)) {
		L = L - w;
		L = max(0, L);
		J = J - 1;
	}

	while (K>0 && z < getPostSigma2(R, eta2[i-1], lambda[i-1], yObs)) {
		R = R+w;
		K = K-1;
	}


	# Shrinkage to obtain a sample
	u = runif(1);
	newParam = L + u*(R-L);
	while (z > getPostSigma2(newParam, eta2[i-1], lambda[i-1], yObs)) {
		if (newParam < currentParam) {
			L = newParam;
		} else {
			R = newParam;
		}
		u = runif(1);
		newParam = L + u*(R-L);
	}

	sigma2[i] = newParam;


# Sampling eta2

	w = 1; m = 20;
	 
	z = getPostEta2(sigma2[i], eta2[i-1], lambda[i-1], yObs) - rexp(1);
	currentParam = eta2[i-1];

	# Stepping out to obtain the [L, R] range
	u = runif(1);
	L = currentParam - w*u;
	R = L + w;
	v = runif(1);
	J = floor(m*v);
	K = (m-1) - J;

	L = max(0, L);
	while (J>0 && z < getPostEta2(sigma2[i], L, lambda[i-1], yObs)) {
		L = L - w;
		L = max(0, L);
		J = J - 1;
	}

	while (K>0 && z < getPostEta2(sigma2[i], R, lambda[i-1], yObs)) {
		R = R+w;
		K = K-1;
	}


	# Shrinkage to obtain a sample
	u = runif(1);
	newParam = L + u*(R-L);
	while (z > getPostEta2(sigma2[i], newParam, lambda[i-1], yObs)) {
		if (newParam < currentParam) {
			L = newParam;
		} else {
			R = newParam;
		}
		u = runif(1);
		newParam = L + u*(R-L);
	}

	eta2[i] = newParam;
	

# Sampling lambda

	w = 0.5; m = 20;

	z = getPostLambda(sigma2[i], eta2[i], lambda[i-1], yObs) - rexp(1);
	currentParam = lambda[i-1];

	# Stepping out to obtain the [L, R] range
	u = runif(1);
	L = currentParam - w*u;
	R = L + w;
	v = runif(1);
	J = floor(m*v);
	K = (m-1) - J;

	L = max(0, L);
	R = min(1, R);

	while (J>0 && L>0 && z < getPostLambda(sigma2[i], eta2[i], L, yObs)) {
		L = L - w;
		L = max(0, L);
		J = J - 1;
	}

	while (K>0 && R<1 && z < getPostLambda(sigma2[i], eta2[i], R, yObs)) {
		R = R+w;
		R = min(1, R);    
		K = K-1;
	}

	# Shrinkage to obtain a sample
	u = runif(1);
	newParam = L + u*(R-L);
	while (z > getPostLambda(sigma2[i], eta2[i], newParam, yObs)) {
    		if (newParam < currentParam) {
			L = newParam;
		} else {
			R = newParam;
		}
    
		u = runif(1);
		newParam = L + u*(R-L);
	}

	lambda[i] = newParam
	
}
return(list(sigma2 = sigma2, eta2 = eta2, lambda = lambda));
}





# The following three functions evaluate the psoterior probability for the most recent parameter values. Note that when updating mu, the prior probability of other parameters are treated as constant. The same is true for eta and lambda.
getPostSigma2 = function(sigma2, eta2, lambda, yObs) {
	prior = log(dgamma(sigma2, alpha1, beta1));
	likelihood = sum(log(lambda*dnorm(yObs, 0, sqrt(sigma2))+ (1-lambda)*dnorm(yObs, 0, sqrt(sigma2 + eta2))));
	return(prior+likelihood);
}

getPostEta2 = function(sigma2, eta2, lambda, yObs) {
	prior = log(dgamma(eta2, alpha1, beta1));
	likelihood = sum(log(lambda*dnorm(yObs, 0, sqrt(sigma2))+ (1-lambda)*dnorm(yObs, 0, sqrt(sigma2 + eta2))));
	return(prior+likelihood);
}


getPostLambda = function(sigma2, eta2, lambda, yObs) {
	prior = log(dbeta(lambda, alpha2, beta2));
	likelihood = sum(log(lambda*dnorm(yObs, 0, sqrt(sigma2))+ (1-lambda)*dnorm(yObs, 0, sqrt(sigma2 + eta2))));
	return(prior+likelihood);
}

# Running the program
res = assignment2Slice(10000);

lambda = mcmc(res$lambda)
sigma2 = mcmc(res$sigma2)
eta2 = mcmc(res$eta2)

plot(lambda)
plot(sigma2)
plot(eta2)

# After convergence diagnosis and burning in some initial states:

summary(lambda)
summary(sigma2)
summary(eta2)

