# This program uses the Metropolis algorithm for the binomial model with Beta prior
# Note that this is a conjugate situation so we can compare our results
# with the result based on the closed form.

# Global Parameters
# These are the parameters of the prior Beta(1, 1).
alpha = 1; beta = 1;

# This the observed data
yObs = 39;
n = 100;

binoPostMetropolis = function() {

	nIterations = 20000
	burnIn = 2000;

	# This is our initial point
	theta = numeric(0);
	theta[1] = 0.5;


	for(i in 2:nIterations) {
		
		# Proposing a new point from Uniform(0, 1) distribution
		thetaPrime = runif(1);

		# Evaluating the posterior probability of the new point
		postThetaPrime = getPosterior(thetaPrime);

		# Evaluating the posterior probability of the current point
		postThetaCurrent = getPosterior(theta[i-1]);

		# Evaluating the acceptance probability
		acceptanceRate = min(1, postThetaPrime/postThetaCurrent);

		# print(thetaPrime, postThetaPrime, 

		u = runif(1);
		if (u<= acceptanceRate) {
			theta[i] = thetaPrime;
		} else {
			theta[i] = theta[i-1];
		}

	}
	
	# Comparing the result with what we obtain using the closed form of the posterior distribution 
	# based on the conjugacy.

	hist(theta[(burnIn + 1):length(theta)], breaks = seq(0, 1, by = 0.02), freq=FALSE, xlab=expression(theta), main='');
	# 8 is a scale factor to make densities line up
	lines(seq(0, 1, by = 0.01), dbeta(seq(0, 1, by = 0.01), shape1 = 40, shape2 = 62), lty=2, lwd=2)

	return(theta);
}

# The following function evaluates the psoterior probability for a given
# theta

getPosterior = function(theta) {

	# We use the following to get the posterior distribution. 
	
	# FOR LARGER SAMPLES, IT IS BETTER TO USE THE LOG TRANSFORMATION OF THESE DISTRIBUTION AND MODIFY THE ACCEPTANCE PROBABILITY ACCORDINGLY.

	prior = theta^(alpha-1)*(1-theta)^(beta-1);
	likelihood =  theta^(yObs) * (1-theta)^(n-yObs);
	posterior = prior*likelihood;
	return(posterior);

	# Alternatively, we could simplify the form of the posterior distribution posterior = (theta^(39))*(1-theta)^(61);
}

theta = binoPostMetropolis();
