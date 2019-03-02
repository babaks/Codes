# This is an example of a good importance sampling. Here we want to use t_3
# for sampling from N(0, 1) so we can approximate the expectation of x^2 
# (with respect to the normal distribution)

	
goodImportanceSampling <- function() {

	# This part plots N(0, 1) and t(3)
	n = 10000;
	dt3 = function(x) dt(x, df = 3);
	plot(dnorm, from = -10, to = 10);
	curve(dt3, from = -10, to = 10, add = TRUE, lty = 2);
	legend(x = 0.4, y = -10, lty = c(2, 1), legend = c("t(3)", "Normal"));

	w = numeric(0);
	samp = numeric(0);

	for(i in 1:n) {
		# We sample from the t distribution with 3 df
		x = rt(1, df = 3);

		# Calculate the unnormalized density of N(0, 1)
		fx = dnorm(x);

		# Calculate the unnormalized density of t_3
		gx = dt(x, df = 3);

		# Calculate the weight
		w[i] = fx / gx;
		samp[i] = x;
	}

	# We use the Monte Carlo formula to estimate E(x^2) which we expect to be very close to its true value 1.
	Ex2 = sum(samp^2 * w / sum(w));
	return(Ex2);
}
