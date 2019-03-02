# This function simulate data for housing prices in 5 states, and uses the
# Gibbs sampler for the hierarchical model described in the session 9.
# The outputs are the true values of the parameters, and the posterior samples 
# for mu (the mean parameter in each state)
# sigma2 (overal variance of housing prices for the US, i.e., not state
# specific), and mu0 (the overall mean of the housing prices in the US).
 

modelHousePrice = function() {

# This part creates the data which includes 5 states each with 100 samples.  
 
# The number of states
J = 5; 
 
# Set the overall average of housing prices (log scale) in the US to 13.37 (i.e.,
# around 640,000$). This would the true value of the parameter. We expect
# the posterior expectation of mu0 based on our MCMC to be close to
# this true value. Of course, since the sample size for mu0 is only 5 (we
# observe only 5 states), the estimates obtained either based on maximum
# likelihood or posterior expectation might differ a little form the true
# value.
mu0_true = 13.37; 
 
# Number of houses samples in each state
nj = 100; 
 
# Using the overall average, and assuming  mu_j ~ N(mu0, 4^2), I sample the
# average housing prices for 5 states. This would the true value of the parameter mu. We expect
# the posterior expectation of mu based on our MCMC to be very close to
# this true value especially since the sample size is large (100).
mu_true = rnorm(J, mu0_true, 4);
# This is the true variance of housing prices which is common between all
# states. We expect the posterior expectation of sigma2 based on our MCMC to be very close to
# this true value especially since the sample size is large (500).
sigma2_true = 5.35;
 
# Withing each state "j", I sample 100 values from N(mu_j, 5.35^2).

state = rep(0, J*nj);
y = rep(0, J*nj);
for(j in 1:J) {
	state[seq((j-1)*nj +1, j*nj, 1)] = j;
	y[seq((j-1)* nj +1, j* nj, 1)] = rnorm(nj, mu_true[j], sqrt(sigma2_true));
	}


# Now that we generated the data, we use MCMC to sample from the posterior
# distributions.
 
# I assume y_ij ~ N(mu_j, sigma2), that is, while the average housing price
# changes from one state to another, the variance remain the same and equal
# to the overal variance of housing prices in the US.
 
 
# The hierarchical model is defined as follows:
 
# y_ij ~ N(mu_j, sigma2)
 
# mu_j ~ N(mu0, ta02)
# sigma2 ~ Inv-chi2(nu0, sigma02)
 
# mu0 ~ N(M, V)
 
 
# Number of MCMC iterations
nIterations = 10000;
 
# mu_j ~ N(0, 25^2), this is the mean of housing prices in each state. As
# you can see, mu0 is a vector (i.e., it's changing), whereas tau02 is
# fixed. mu0 is a hyperparameter here.
mu0 = rep(0, nIterations);
tau02 = 25^2;
 
# mu0 ~ N(0, 50^2), the parameters of the normal distribution for 
M = 0;
V = 50^2;
 
 
# sigma2 ~ Inv-chi2(nu0, sigma02), note that nu0 and and sigma02 are fixed,
# they are not hyperparameters.
nu0 = 1;
sigma02 = 0.5;
 
# These are the parameters of y_ij ~ N(mu_j, sigma2). Both of them are
# random variables, however, there is only 1 sigma2, i.e., it's common between all
# states, whereas there are 5 different mu's, one for each state.
mu = matrix(rep(0, nIterations*J), nrow=nIterations, ncol = J)
sigma2 = matrix(rep(1, nIterations), nrow=nIterations, ncol = 1)
 

for(i in 2:nIterations) {
    # yCent is a variable that holds (y_ij - mu_j), I need this to
    # calculate nu (which is needed to get the posterior distribution of
    # sigma2). 
    yCent = rep(0, nj*j);
    
    # Here, mu_0 and simga2 are fixed. I sample mu_j's one at a time. 
    for(j in 1:J){
        
    # Housing prices for state j.    
    yObs = y[state==j];
    # samples size in that state.
    n = length(yObs);
    
    # This is (conditionally) conjugate situation, so I can directly sample
    # from the posterior distribution mu_j ~ N(mu_n, sigma2_n). 
    
    # Calculating mu_n and sigma2_n: these formulas are the same as what we saw for simple
    # normal model with conjugate prior and known variance.
    mu_n = (mu0[i-1]/tau02 + sum(yObs)/sigma2[i-1]) / (1/tau02 + n/sigma2[i-1]);
    sigma2_n = 1/(1/tau02 + n/sigma2[i-1]);
    
    # Sampling a new value for mu_j at iteration i.
    mu[i, j] = rnorm(1, mu_n, sqrt(sigma2_n));
    yCent[seq((j-1)*nj +1, j*nj, 1)] =  yObs - mu[i, j];
    }

    # Since sigma2 is common for all states, I use all the housing prices
    # regardless of their states to obtain its posterior distribution. Note
    # that when I was subtracting the mean from y to calculate nu, I
    # subtracted mu_j from y_ij (i.e., for each state separately) because 
    # these y's now have different means.
    
    # Given mu_j, which we sampled in the last step, this becomes a (conditionally)
    # conjugate situation and similar to normal models with known mean.
    n = length(y);    

	nu_n = nu0+n;
	nu = sum(yCent^2)/n;
	sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
	
    # I sample a new value for sigma2 from Inv-chi2(mu_n, sigma02_n), to do
    # this, I sample from z ~ chi2(nu_n) and then my sample from Inv-chi2
    # would be nu_n*sigma02_n/z.

	
	z = rchisq(1, nu_n);
	sigma2[i] = nu_n*sigma02_n/z;


    # Now, using the new values for mu_j, I go one level up in the
    # hierarchy, and sample a new mu_0. Note that this is the same as what
    # we used for a normal model with conjugate prior for mean and
    # known variance. The only difference here now is that our data are not
    # the observed housing prices, they are the new samples for mu_j.
    
    # Therefore, our sample size is 5 since we have 5 states.
    n = J;
    
    # I sample from the posterior distribution mu0 ~ N(mu_n, sigma2_n)
    mu_n = (M/V + sum(mu[i, ])/tau02) / (1/V + n/tau02);
    sigma2_n = 1/(1/V + n/sigma02);
    mu0[i] = rnorm(1, mu_n, sqrt(sigma2_n));


	}
	return(list(mu0_true = mu0_true, mu_true = mu_true, sigma2_true = sigma2_true, mu = mu, sigma2 = sigma2, mu0 = mu0));
}

param = modelHousePrice();