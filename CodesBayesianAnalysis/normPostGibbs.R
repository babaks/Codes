# This program uses the Gibbs sampler for simulating from the posterior
# distribution of a normal model with unknown mean and variance. The mean
# and variance are given normal and Inv-chi2 priors. 

par(ask = TRUE);
normPostGibbs = function() {

	nIterations = 200;

	mu0 = 65;
	tau02 = 100^2;

	nu0 = 1;
	sigma02 = 1;

	yObs = c(72, 75, 70);
	n = length(yObs);

	mu = 50;
	sigma2 = 1;

	par(mfcol = c(1, 2));
	plot(mu, mfg = c(1, 1), xlab = "iteration", ylab =  "mu", type = 'l');
	plot(sqrt(sigma2), mfg = c(1, 2), xlab = "iteration", ylab =  "sigma", type = 'l');

	for(i in 2:nIterations) {
	# The conditional distribution of mu given sigma is normal(mu_n,
	# sigma_n) when we use a normal(mu0, tau02) prior for mu.
	mu_n = (mu0/tau02 + sum(yObs)/sigma2[i-1]) / (1/tau02 + n/sigma2[i-1]);
	sigma_n = sqrt(1/(1/tau02 + n/sigma2[i-1]));
	mu[i] = rnorm(1, mean = mu_n, sd = sigma_n);

	# The conditional distribution of variance sigma2 is Inv-chi2(nu_n, sigma02_n)
	# We can first sample z form chi2(nu_n) distribution then use nu_n*sigma02_n/z
	# as a sample from the scaled Inv-Chi2 distribution.  
	nu_n = nu0+n;
	nu = (sum( (yObs - mu[i])^2 ))/n;
	sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
	z = rchisq(1, nu_n);
	sigma2[i] = nu_n*sigma02_n/z;

		par(mfcol = c(1, 2));
		plot(mu, mfg = c(1, 1), xlab = "iteration", ylab =  "mu", type = 'l');
		plot(sqrt(sigma2), mfg = c(1, 2), xlab = "iteration", ylab =  "sigma", type = 'l');

	}
	return(list(mu = mu, sigma2 = sigma2));
}

param = normPostGibbs();