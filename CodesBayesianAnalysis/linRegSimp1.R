# This is a Bayesian linear regression model with informative priors The outputs are MCMC samples for the regression parameters "beta" and sigma2, the variance of (y|x, beta). We are using a simple Gibbs sampler here. The program could be written a little more efficiently, but we are trying to keep a balance between efficiency of the code and its readability.  

# Number of MCMC iterations
nIter = 2000;
n = 2000 
p = 1

# Reading in the data
x = c(rep(0, n/2), rep(1, n/2))
y = -1+2*x + rnorm(n, 0, sqrt(1.5))

# We choose mother's IQ and whether she is graduated from high school as predictors 
# The output is kid's test score

 
# Here, we are standardizing the predictors by subtracting their mean and deviding by their standard deviation.
# x = scale(x); 

# Adding a column of 1's for the intercept
x = cbind(rep(1, n), x); 

# d is the number of regression parameters including the intercept.
d = dim(x)[2]
 
# The parameteres of Inv-chi2 prior for sigma2
nu0 = 1;
sigma02 = 0.5;
 
# The parameters of normal prior N(mu0, tau02) for beta's.
mu0 = rep(0, d);
tau02 = rep(100, d);
  
# The following two matrices hold the posterior samples for sigma2 and beta.
sigma2 = rep(1, nIter);
beta = matrix(0, nIter, d);
 
# The prior distribution for beta in the form of a matrix. 
Lambda0 = diag(tau02);
invLambda0 = diag(1/tau02);
 
 
# It is always better to use the Cholosky decomposition for getting the iverse of matrices and sampling from multivariate normal. Here, we first get, L1, the Cholesky decomposition of (x'x). To get, the inv(x'x), we use chol2inv.
 
X = t(x)%*%x;
L1 = chol(X);
invX = chol2inv(L1);
 
# This give the mean of the posterior distribution.
# betaHat = invX%*%t(x)%*%y;
 
 
# These are y*, x* as described in the course notes.
yStar = c(y, mu0);
xStar = rbind(x, diag(d));
  
for(i in 2:nIter){
    
    # A contains the diagonal elements of SigmaStar^(-1/2)
    A = 1/sqrt(c(sigma2[i-1]*rep(1, n), tau02));
    # newX is SigmaStar^(-1/2)*xStar
    newX = matrix(rep(A, d), nrow=(n+d), ncol = d)*xStar
        
    # newY is SigmaStar^(-1/2)*yStar
    newY = A*yStar;
    
    # Now we can get mu_n and Lambda_n as described in the course note based
    # on newX and newY
    L2 = chol(t(newX)%*%newX);
    invL2 = solve(L2);
    Lambda_n = invL2%*%t(invL2);
    mu_n = Lambda_n%*%t(newX)%*%newY;
    
    # We can now sample from the posterior distribution given sigma^2. As
    # discussed in the class, we use the Cholesky decomposition for this.
    L3 = chol(Lambda_n);
    u = rnorm(d, 0, 1);
    beta[i, ] = u%*%L3 + t(mu_n);
 
	# Now, given the current beta's, we sample from the posterior distribution of 	sigma2, which is Inv-chi2(nu_n, sigma02_n). 
    eta = x%*%beta[i, ];
    eps = y - eta;
    nu = sum(eps^2)/(n);
    nu_n = nu0 + n;
    sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
    
    # I sample a new value for sigma2 from Inv-chi2(nu_n, sigma02_n), to do this, I 	sample from z ~ chi2(nu_n) and then my sample from Inv-chi2 would be 	nu_n*sigma02_n/z.
    z = rchisq(1, nu_n);
    sigma2[i] = nu_n*sigma02_n/z;

   
    }

