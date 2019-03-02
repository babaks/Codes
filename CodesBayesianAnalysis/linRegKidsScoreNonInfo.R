# This is a Bayesian linear regression model with noninformative priors. The outputs are MCMC samples for the regression parameters "beta" and sigma2, the variance of (y|x, beta). We are using a simple Gibbs sampler here.


linRegKidsScoreNonInfo <- function(){

# Number of MCMC iterations 
nIterations = 10000;
 
# Reading in the data
data = as.matrix(read.table("kidiq.txt", sep = ",", header=TRUE));

# We choose mother's IQ and whether she is graduated from high school as predictors 
# The output is kid's test score
x = data[, 2:3];
y = data[, 1];


n = dim(x)[1]
p = dim(x)[2] 
 
# Here, we are standardizing the predictors by subtracting the mean and deviding by their standard deviation.
x = scale(x); 

#Adding a column of 1's for the intercept
x = cbind(rep(1, n), x);
 
 
# The following two matrices hold the posterior samples for sigma2 and beta.
sigma2 = rep(1, nIterations);
beta = matrix(rep(0, nIterations*(p+1)), nrow = (p+1), ncol=nIterations);

 
# It is always better to use the Cholosky decomposition to invert matrices and sampling from multivariate normal. Here, we first get, L1, the Cholesky decomposition of (x'x). To get, the inv(x'x), we use chol2inv function.

X = t(x)%*%x;
L1 = chol(X);
invX = chol2inv(L1);
 
# This give the mean of the posterior distribution.
betaHat = invX%*%t(x)%*%y;
 
# Since later on we need the sample from N(betaHat, invX*sigma2), we should find the Cholesky decomposition of invX*sigma2. Instead of doing this in a loop, we can get the Cholesky of invX (we call it L2) and then at each iteration, we multiply this by sqrt(sigma2).
L2 = chol(invX);
 
for(i in 2:nIterations){
    
    L3 = L2*sqrt(sigma2[i-1]);       
    


	# We can now sample from the posterior distribution of betea 	given sigma^2. As discussed in the class, we use the Cholesky 	decomposition for this.
    
    u = rnorm(p+1, 0, 1);
    beta[, i] = u%*%L3 + t(betaHat);
    
 
    # Now, given the current beta's, we sample from the posterior
    # distribution of sigma2, which is Inv-chi2(n-p-1, s2).
    eta = x%*%beta[, i];
    eps = y - eta;
    s2 = sum(eps^2)/(n-p-1);

    # I sample a new value for sigma2 from Inv-chi2(n-p-1, s2). To 	do this, I sample from z ~ chi2(n-p-1) and then my sample from 	Inv-chi2 would be (n-p-1)*s2/z.    
    z = rchisq(1, n-p-1);
    sigma2[i] = (n-p-1)*s2/z;
    
    }
	
	return(list(beta = beta, sigma2 = sigma2))

}


# Running the program
res<-linRegKidsScoreNonInfo()



