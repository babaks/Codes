# note: packages must be installed from the CRAN website before they can be used in R;
# this can be done directly from the R GUI (Packages --> Install Packages in the Windows version)

par(ask = TRUE);
# load the necesarry packages
# need this package for multivariate normal
library(mvtnorm);

# This program uses the Metropolis algorithm for the multivariate normal 
# model with unknown mean and known covariance matrix. We use a
# multivariate normal prior for the mean.
# Note that this is a conjugate situation so we can compare our results
# with the result based on the closed form.
# The difference between this program and normPostMetroplis1.m is that this
# oen updates parameters one at a time. 

normPostMetropolis2 <- function() {
	# Number of iterations and number of samples we discard (burn in).
	nIterations = 1000;
	burnIn = 100;

	# These are the parameters of the prior
	mu0 = c(0, 0);
	lambda0 <- diag(c(10, 10));

	yObs <- rbind(c(-1.2, 2.3), c(-0.5, 0.7), c(-2.1, -1));
	n = dim(yObs)[1];
	d = dim(yObs)[2];

	sigma = cov(yObs);

	theta = matrix(0, nrow = nIterations, ncol = 2);

	mu_n = (mu0 %*% solve(lambda0) + n*apply(yObs, 2, mean) %*% solve(sigma)) %*% solve(solve(lambda0)+n*solve(sigma));
	lambda_n = solve(solve(lambda0)+n*solve(sigma));

	x = seq(-4, 4, by = 0.2);
	y = seq(-4, 4, by = 0.2);
	len = length(x);
	thetaConjugate = matrix(0, nrow = len, ncol = len);
	for(i in 1:len) {
		for(j in 1:len) {
			thetaConjugate[i, j] = dmvnorm(x = c(x[i], y[j]), mean = mu_n, sigma = lambda_n);
		}
	}
	contour(x = x, y = y, z = thetaConjugate, main = "The contour plot shows the conjugate posterior distribution", 
		xlab = expression(mu[1]), ylab = expression(mu[2]));
	points(theta[1, 1], theta[1, 2]);

	for(i in 2:nIterations) {
		thetaPrime = theta[i-1,];
		for(j in 1:d) {
			# Proposing one dimension of a new point from a Normal distribution
			thetaPrime[j] = rnorm(n = 1, mean = theta[i - 1,j], sd = .5);
 
			# Evaluating the posterior probability of the new point
    			postThetaPrime = getPosterior(thetaPrime, yObs, n, sigma, mu0, lambda0);
    
    			# Evaluating the posterior probability of the current point
    			postThetaCurrent = getPosterior(theta[i-1, ], yObs, n, sigma, mu0, lambda0);
    
    			# Evaluating the acceptance probability
    			acceptanceRate = min(1, postThetaPrime/postThetaCurrent);
    
    			# Deciding whether to accept the new point or stay where we are
    			u = runif(1);
    			if (u<= acceptanceRate) {
        			theta[i, j] = thetaPrime[j];
    			} else {
        			theta[i, j] = theta[i - 1, j];
    			}

			if (i <= 20) {
				contour(x = x, y = y, z = thetaConjugate, main = "The contour plot shows the conjugate posterior distribution", 
					xlab = expression(mu[1]), ylab = expression(mu[2]));
				if(j == 1) {
					lines(rep(theta[1:i,1], times = c(1, rep(2, i - 2), 1)), 
						rep(theta[1:(i-1), 2], times = rep(2, i - 1)));
				} else {
					lines(rep(theta[1:i,1], times = c(1, rep(2, i - 1))), 
						rep(theta[1:i, 2], times = c(rep(2, i - 1), 1)));
				}
			}
		}
	}

	contour(x = x, y = y, z = thetaConjugate, main = "The contour plot shows the conjugate posterior distribution", 
		xlab=expression(mu[1]), ylab = expression(mu[2]));
	points(theta[burnIn:nIterations, 1], theta[burnIn:nIterations, 2]);

	plot(theta[,1], col="blue", type='l',
		xlab = "Iteration", ylab = expression(mu), ylim = c(-2.5, 2.5));
	lines(theta[,2], col = "red");
	legend(legend = c(expression(mu[1]), expression(mu[2])), col = c("blue", "red"), x = 0, y = 2.5, lty = c(1, 1));

	return(theta);
}

getPosterior = function(theta, yObs, n, sigma, mu0, lambda0) {
	prior = dmvnorm(theta, mu0, lambda0);
	likelihood = prod(dmvnorm(yObs, theta, sigma));
	posterior = prior * likelihood;
	return(posterior);
}





par = normPostMetropolis2();