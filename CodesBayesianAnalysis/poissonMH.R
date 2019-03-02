# options(error = recover);

# This program uses the Metropolis-Hastings (MH) algorithm for the Poisson 
# model with Gamma prior Note that this is a conjugate situation so we can 
# compare our results with the result based on the closed form.

poissonPostMH = function() {
	# Number of iterations and number of samples we discard (burn in).
	nIterations = 20000;
	burnIn = 2000;

	# These are the parameters of the prior
	alpha = 1.4; beta = 10;

	yObs = c(0, 1);
	n = length(yObs);

	theta = 0.5;

	for(i in 2:nIterations) {
		# Proposing a new point from the Uniform(0, theta+1) distribution distribution
		thetaPrime = runif(1, 0, theta[i-1]+1);

		# Evaluating the posterior probability of the new point
   		postThetaPrime = getPosterior(thetaPrime, yObs, n, alpha, beta);
    
    		# Evaluating the posterior probability of the current point
   		postThetaCurrent = getPosterior(theta[i-1], yObs, n, alpha, beta);
		
		# Transition probability from the current theta to theta'
		g.thetaCurrent.thetaPrime = dunif(thetaPrime, 0, theta[i-1]+1);
    
		# Transition probability from theta' to the current theta
		g.thetaPrime.thetaCurrent = dunif(theta[i-1], 0, thetaPrime+1);
    
    		# Evaluating the acceptance probability
 		acceptanceRate = min(1, (postThetaPrime*g.thetaPrime.thetaCurrent)/(postThetaCurrent*g.thetaCurrent.thetaPrime));    


    		# Deciding whether to accept the new point or stay where we are
    		u = runif(1);
    		if (u<= acceptanceRate) {
        		theta[i] = thetaPrime;
    		} else {
        		theta[i] = theta[i - 1];
    		}

	}

	theta0 = seq(0, 2, by = 0.01);
	thetaConjugate = dgamma(theta0, 2.4, 12);
	hist(theta[burnIn:nIterations], xlab = expression(theta), ylab = "density", freq = FALSE, main='')
	lines(theta0, thetaConjugate, lty=2, lwd=2)
	return(theta);
}

getPosterior = function(theta, yObs, n, alpha, beta) {
	prior = theta^(alpha-1)*exp(-beta*theta);;
	likelihood = theta^(sum(yObs))*exp(-n*theta);
	posterior = prior * likelihood;
	return(posterior);
}





theta = poissonPostMH();