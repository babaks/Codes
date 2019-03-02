# This is the computer program for the rats tumor model we discussed in the
# class. Here, we use a hierarchical model where we use one parameter
# theta_j for the probability of developing tumor in each rats group.
# Moreover, we assume all theta's share a common Beta(alpha, beta)
# distribution. For simplisity, we fix alpha, and regard beta as a hyperparameter. 
# The hyperprior distribution for beta is Gamma(gamAlpha, gamBeta).
# The outputs of this program are posterior samples for each theta, and the
# posterior samples for beta. The outputs are param$beta and param$theta.
 
nIterations = 10000;
 
# Reading in the data.
rats = read.table('ratTumor.txt', sep = ",");
y = rats[, 1];
n = rats[, 2];
nRats = length(y);
J = 71;
 
# These are the maximum likelihood estimates for theta's.
thetaObs = y/n;
 
# Setting up the parameters of the priors.
gamAlpha = 2;
gamBeta = 1;
alpha = 2;

# These two matrices hold posterior samples 
theta = matrix(rep(0, nIterations*J), nrow=nIterations, ncol = J); 
beta = rep(1, nIterations);

modelRats = function(){ 
 
 
for (i in 2:nIterations){
    # Given the current beta, the posterior distribution of each theta_j
    # has a closed form (Beta distribution) since the prior is conjugate. 
    # So we can use the Gibbs sample to update all the 71 theta's. As you 
    # can see, for each of them, we use their own specific y_j and n_j. 
    for (j in 1:nRats){
        theta[i, j] = rbeta(1, alpha+y[j], beta[i-1]+n[j]-y[j]);
    }
browser()    
    # Now we move one level higher and take theta's as observations and
    # update beta using a Metropolis step (note, we could use conjugacy and
    # update beta using the Gibbs sampler, but I haven't talked about this
    # in class, so I use Metropolis).
    # We use Uniform(max(0, beta-0.5), beta+0.5). As you can see, whenever 
    # theta-0.5 becomes negative, I use 0 as the lower bound since beta has
    # to be positive and it is not efficient to propose negative values
    # which would be rejected. This is still symmetric proposal since
    # uniform distribution has a constant density.
    lower = max(0, beta[i-1]-1);
	upper = beta[i-1]+1;
	
    gCurrent = 1/(upper-lower);

    betaPrime = runif(1, lower, upper);
     
    lower = max(0, betaPrime-1);
    upper = betaPrime+1;
    
    gPrime = 1/(upper-lower);

    
    # This obtains the log-posterior probabilities for the proposed sample
    # and the current sample.
    logPostPrime = getPost(betaPrime, theta[i, ]);
    logPostCurrent = getPost(beta[i-1], theta[i, ]);
    
    # Note that I am using the log transformation of posterior
    # distribution, so I am modifying the accceptance probability.
    acceptanceProb = min(1, exp(logPostPrime + log(gPrime) - logPostCurrent - log(gCurrent)));
 
    # Accepting new proposed values with probability acceptanceProb.
    u = runif(1);
    if (u<= acceptanceProb){
        beta[i] = betaPrime;
    } else {
        beta[i] = beta[i-1];
    }
 browser()
}
    # This gives the posterior expectation of the mean for Beta distribution, which can be 
    # interpreted as the overall mean (over all rat groups) of developing tumor.  
    overallMean = mean(alpha/(alpha+beta[seq(1000, nIterations, 1)]));
	
    return(list(theta=theta, beta=beta, overallMean = overallMean));
}
  
# This is the part we calculate the log-posterior probability for a given
# beta and theta. The hyperprior for beta is gamma and the p(theta|alpha,
# beta)=Beta(alpha, beta).
getPost = function(beta, theta){
logPrior = log(dgamma(beta, gamAlpha, gamBeta));
logLike = sum(log(dbeta(theta, alpha, beta)));
logPost = logPrior + logLike;
return(logPost)
}


param = modelRats();

