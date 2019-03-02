# This program performs rejection sampling from Beta(3, 10) using a Uniform(0, 1) and c=4.

rejectionSampling <- function() {

	par(ask = TRUE);


	# Number of samples
	n = 10000;

	# Ploting the actual distribution for comparision and to see that 4*Uniform(0, 1) has the blanket property.

	x0 = seq(0, 1, by = 0.01);
	y1 = dbeta(x0, shape1 = 3, shape2 = 10);
	plot(x0, y1, ylim = c(0, 5), type = 'l');
	y2 = rep(4, length(x0));
	lines(x0, y2, lty = 2);
	legend(x = 0, y = 5, lty = c(1, 2), legend = c("Beta(3, 10) distribution", "Uniform(0, 1) multiplied by c = 4"));

	samp = numeric(0);

	for(i in 1:n) {
		# We first sample from Uniform(0, 1)
		x = runif(1);
		
		# We then evaluate f*(x). I am of course using f(x) here for illustration purpose. In general, we can drop the constant.
		f.x = dbeta(x, 3, 10)
		# We sample again from Uniform(0, cg*(x)) but this time to decide whether we should accept or reject the sample. Here, we set c=4.
		u = runif(1, 0, 4);
		if (u<=f.x) {
		samp = c(samp, x);
		}
	}
	
	plot(x0, y1, ylim = c(0, 5), type = 'l');
	lines(x0, y2, lty = 2);
	f = hist(samp, breaks = seq(0, 1, by = .02), plot = FALSE);
	lines(f$breaks[-length(f$breaks)], 50 * f$counts /length(samp), type = 's');
	legend(x = 0, y = 5, lty = c(1, 2), legend = c("Beta(3, 10) distribution", "Uniform(0, 1) multiplied by c = 4"));	
	
}
