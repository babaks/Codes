# In this example, we consider a normal model with known variance. We use normal prior for the mean and use MCMC for sampling from the posterior distribution. We compare this distribution with the colsed form based on the conjugacy of prior.

# The inputs are the number of iterations (nIterations) and the step size (delta), which here is the varaince of the proposal distribution.

normMetropolis = function(nIterations = 1000, delta = 1) {
	
	# These are the parameters of the prior normal distribution for the mean
	mu0 = 65; tau0 = 3

	# This is the observed data
	yObs = c(72, 75, 70)
	
	n = length(yObs);
	sigma = 4 # standard deviation of the data, which is assumed to be known
	
	# theta is a vector that stors the posterior samples
	theta = rep(NA, nIterations)
	
	# We use the sample mean as our starting point
	theta[1] = 60
	
	for(i in 2:nIterations) {
		
		# Proposing a new point from the N(theta, delta) distribution distribution
		thetaPrime = rnorm(1, theta[i-1], delta)
		
		# Evaluating the posterior probability of the new point
   		postThetaPrime = getPosterior(thetaPrime, yObs, sigma, mu0, tau0);
    
    	# Evaluating the posterior probability of the current point
   		postThetaCurrent = getPosterior(theta[i-1], yObs, sigma, mu0, tau0);
		
		# Calculating the acceptance probability 
		acceptanceRate = min(1, postThetaPrime/postThetaCurrent);

		# Deciding whether to accept the proposed value or not
		u = runif(1);
		if (u<= acceptanceRate) {
			theta[i] = thetaPrime;
		} else {
			theta[i] = theta[i-1];
		}


	}

	# Returning posterior samples
	return(theta);
}


# This function calculates the posterior density up to a constant
getPosterior = function(theta, yObs, sigma, mu0, tau0) {
	prior = exp(-(theta - mu0)^2/(2*tau0^2))
	likelihood = exp(- sum((yObs-theta)^2)/(2*sigma^2))
	posterior = prior * likelihood
	return(posterior)
}



# Running the program
nIterations = 10000
delta = 1
theta = normMetropolis(nIterations, delta)

par(mfrow=c(2, 1))
plot(theta, xlab='Iteration', ylab=expression(theta), type='l')

# Discarding the first 5000 samples
burnIn = 5000

# Plotting the histogram of posterior samples
hist(theta[burnIn:nIterations], xlab = expression(theta), ylab = "density", freq = FALSE, main='')

# Comparing it with the posterior distribution from the closed form based on the conjugacy of prior
theta0 = seq(50, 90, by = 0.01);
thetaConjugate = dnorm(theta0, 69.6, sqrt(3.4));
lines(theta0, thetaConjugate, lty=2, lwd=2)