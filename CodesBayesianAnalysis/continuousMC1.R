# This is a continuous space Markov chain with N(x/2, sqrt(3/4)) transition
# distribution. R considers the second parameter as standard
# deviation not variance.

options(error = recover);

par(ask = TRUE);
samp = numeric(0);

continuousMC1 = function() {

	# Number of samples from the Markov chain
	n = 1000;

	X = 0.5;
	for(i in 1:n) {
		# when we are at point X, we sample a new point from N(X/2, sqrt(3/4)).
		X = rnorm(1, mean = X/2, sd = sqrt(3/4));
		samp[i] = X;

		# This part plots the samples. and shows where the transition
		# distributon is at any given time.
		plot(X, 0, xlim = c(-4, 4), ylim = c(-0.05, 0.6), pch = 2, main = "Sampling from a Markov chain with N(X/2, 3/4) transition distribution");
		lines( seq(-4, 4, by=.1), dnorm(seq(-4, 4, by=.1),  mean = X/2, sd = sqrt(3/4)));
		print(i);
		if(i >= 2) {
			points(samp[1:(i - 1)], rep(0, i - 1));
		}
		legend(x = -4, y = .6, pch = c(1, 2), legend = c("Current State", "Previous States"))
	}
	hist(samp[101:n]);
}

continuousMC1();

