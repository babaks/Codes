# This is a Bayesian linear regression model with informative priors The outputs are MCMC samples for the regression parameters "beta" and sigma2, the variance of (y|x, beta). We are using a simple Gibbs sampler here. The program could be written a little more efficiently, but we are trying to keep a balance between efficiency of the code and its readability.  

# Number of MCMC iterations
nIter = 2000;
n = 200 
p = 100

# Reading in the data
x = c(rep(0, n/2), rep(1, n/2))
y = cbind(x, 2*x, 3*x, matrix(rnorm((p-3)*n), n, (p-3))) + matrix(rnorm(n*p), n, p)
c = x+1
n1 = sum(c==1)
n2 = sum(c==2)

# We choose mother's IQ and whether she is graduated from high school as predictors 
# The output is kid's test score

 
# Here, we are standardizing the predictors by subtracting their mean and deviding by their standard deviation.
# x = scale(x); 

# Adding a column of 1's for the intercept
x = cbind(rep(1, n), x); 

# d is the number of regression parameters including the intercept.
d = dim(x)[2]

X = t(x)%*%x;
L.x = chol(X);
invX = chol2inv(L.x);
L.x.inv = chol(invX);
 
# This give the mean of the posterior distribution.
BETA.hat = invX%*%t(x)%*%y;

x.a = matrix(x[, 1:(d-1)], n, d-1, byrow=TRUE)
x.b = matrix(x[, d], n, 1)

X.a = t(x.a)%*%x.a;
L.a = chol(X.a);
invX.a = chol2inv(L.a);
L.a.inv = chol(invX.a);


X.b = t(x.b)%*%x.b;
L.b = chol(X.b);
invX.b = chol2inv(L.b);
L.b.inv = chol(invX.b);


# The following two matrices hold the posterior samples for sigma2 and beta.
alpha = matrix(BETA.hat[1:(d-1), ], nrow=d-1, ncol=p, byrow=TRUE)
beta = matrix(BETA.hat[d, ], 1, p)

sigma2 = colSums((y - x%*%BETA.hat)^2)/(n-1)

post.sigma2 = matrix(NA, nIter, p)
post.alpha = array(NA, c(nIter, d-1, p))
post.beta = matrix(NA, nIter, p)

 
# The parameteres of Inv-chi2 prior for sigma2
nu0 = 1;
sigma02 = 0.5;
 
# The parameters of normal prior N(mu0, tau02) for beta's.
mu0 = 0
tau02 = 100
  
 
# It is always better to use the Cholosky decomposition for getting the iverse of matrices and sampling from multivariate normal. Here, we first get, L1, the Cholesky decomposition of (x'x). To get, the inv(x'x), we use chol2inv.

 
    
for(i in 1:nIter){
print(i)            
z = y - x.b%*%beta

alpha.hat = invX.a%*%t(x.a)%*%z


sigma2.mat = matrix(sigma2, d-1, p, byrow=TRUE)
    
	# We can now sample from the posterior distribution of betea 	given sigma^2. As discussed in the class, we use the Cholesky 	decomposition for this.
    
    u = matrix(rnorm((d-1)*p), p, d-1);
    alpha = (t(u%*%L.a.inv) + alpha.hat/sqrt(sigma2.mat))*sqrt(sigma2.mat)
    
    
	z = y - x.a%*%alpha
	z = z[c==2, ]
	V = 1/(1/tau02 + n2/sigma2)
	mu_n = V*(mu0/tau02 + colSums(z)/sigma2) 
	sigma_n = sqrt(V)
	beta = rnorm(p, mean = mu_n, sd = sigma_n)
	
	BETA = rbind(alpha, beta)
    # Now, given the current beta's, we sample from the posterior
    # distribution of sigma2, which is Inv-chi2(n-p-1, s2).
    eps = y - x%*%BETA
    nu = colSums(eps^2)/(n);
    nu_n = nu0 + n;
    sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
    
    # I sample a new value for sigma2 from Inv-chi2(nu_n, sigma02_n), to do this, I 	sample from z ~ chi2(nu_n) and then my sample from Inv-chi2 would be 	nu_n*sigma02_n/z.
    z = rchisq(p, nu_n);
    sigma2 = nu_n*sigma02_n/z;
	   
	post.alpha[i, , ] = alpha
	post.beta[i, ] = beta
	post.sigma2[i, ] = sigma2
   
    }

