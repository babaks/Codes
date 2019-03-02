# This is an example where the importance sampling will not work properly.
# Here, we are using N(0, 1) to sample from t_3.

badImportanceSampling <- function() {
	par(ask = TRUE);

	# This part plots N(0, 1) and t_3
	n = 1000;
	dt3 = function(x) dt(x, df = 3);
	plot(dnorm, from = -10, to = 10, lty = 2);
	curve(dt3, from = -10, to = 10, add = TRUE);
	legend(x = 0.4, y = -10, lty = c(1, 2), legend = c("t(3)", "Normal"));

	w = numeric(0);
	samp = numeric(0);

	for(i in 1:n) {
		# We sample from N(0, 1)
		x = rnorm(1);

		# Calculate the unnormalized density of t_3
		fx = dt(x, df = 3);

		# Calculate the unnormalized density of N(0, 1)
		gx = dnorm(x);

		# Calculate the weight
		w[i] = fx / gx;
		samp[i] = x;
	}

	# Use the Monte Carlo formula to estimate the expectation of x^2 with
	# respect to t_3 distribution. The answer would be systematically lower
	# than the actual value 3.
	Ex2 = sum(samp^2 * w / sum(w));
	return(Ex2);
}
