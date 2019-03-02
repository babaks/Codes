## The original codes for this analysis were written by Gelman and Hill. I have modified the codes for Stat211. 


# Getting the data and selecting the counties in Minnesota


data <- read.table ("radon.dat", header=T, sep=",")

head(data)

mn <- data$state=="MN"
radon <- data$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- data$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon

# get county index variable
county.name <- as.vector(data$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}



nj <- as.vector (table (county))

# Jittering sample sizes so the plots don't over on top of each other
sample.size.jittered <- nj*exp (runif (J, -.1, .1))




radon = function(nIterations=1000) {

 

# Now that we generated the data, we use MCMC to sample from the posterior
# distributions.
 
# y_ij ~ N(mu_j, sigma2_j),  
# The hierarchical model is defined as follows:
# y_ij ~ N(mu_j, sigma2_j)
# mu_j ~ N(mu0, ta02)
# sigma2_j ~ Inv-chi2(nu0, sigma02)
# mu0 ~ N(M, V)
 
  
mu0 = rep(0, nIterations);
tau02 = rep(1, nIterations);
 
# mu0 ~ N(0, 10^2), the parameters of the normal distribution for 
M = 0;
V = 10^2;
 
 
# sigma2_j ~ Inv-chi2(nu0, sigma02), note that nu0 and and sigma02 are fixed.
# I also assume tau02 ~ Inv-chi2(nu0, sigma02)
nu0 = 1;
sigma02 = 0.5;
 
mu = matrix(rep(0, nIterations*J), nrow=nIterations, ncol = J)
sigma2 = matrix(rep(1, nIterations*J), nrow=nIterations, ncol = J)
 

for(i in 2:nIterations) {

    
    # Here, mu_0 and simga2 are fixed. I sample mu_j's one at a time. 
    for(j in 1:J){
        

    yObs = y[county==j];

    n = length(yObs);
    
    
    # Calculating mu_n and sigma2_n: these formulas are the same as what we saw for simple normal model with conjugate prior and known variance.
    mu_n = (mu0[i-1]/tau02[i-1] + sum(yObs)/sigma2[i-1]) / (1/tau02[i-1] + n/sigma2[i-1]);
    sigma2_n = 1/(1/tau02[i-1] + n/sigma2[i-1]);
    
    # Sampling a new value for mu_j at iteration i.
    mu[i, j] = rnorm(1, mu_n, sqrt(sigma2_n));
  
	nu_n = nu0+n;
	nu = sum((yObs - mu[i, j])^2)/n;
	sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
	
    # I sample a new value for sigma2 from Inv-chi2(mu_n, sigma02_n), to do this, I sample from z ~ chi2(nu_n) and then my sample from Inv-chi2 would be nu_n*sigma02_n/z.
	z = rchisq(1, nu_n);
	sigma2[i, j] = nu_n*sigma02_n/z;
	}
	
	
    # Now, using the new values for mu_j, I go one level up in the
    # hierarchy, and sample a new mu_0 and a new tau2_0. 
    
    n = J;
    
    # I sample from the posterior distribution mu0 ~ N(mu_n, sigma2_n)
    mu_n = (M/V + sum(mu[i, ])/tau02) / (1/V + n/tau02);
    sigma2_n = 1/(1/V + n/sigma02);
    mu0[i] = rnorm(1, mu_n, sqrt(sigma2_n));

	nu_n = nu0+n;
	nu = sum((mu[i, ] - mu0[i])^2)/n;
	sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
	z = rchisq(1, nu_n);
	tau02[i] = nu_n*sigma02_n/z;
	
	
	}
	return(list(mu = mu, sigma2 = sigma2, mu0 = mu0));
}

nIter = 5000
param = radon(nIter);





########

# Complete pooling
ybarbar = mean(y)

# county means, variance, and standard deviation
cty.mns = tapply(y,county,mean)
cty.vars = tapply(y,county,var)
cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size)
cty.sds.sep = sqrt(tapply(y,county,var)/sample.size)


# no predictors
## no-pooling vs. pooling
par(mfrow=c(1,2))
par (mar=c(5,5,4,2)+.1)
plot (sample.size.jittered, cty.mns, cex.lab=1.1, cex.axis=1.1,
      xlab="sample size in county j",
      ylab="avgerage log radon in county j",
      pch=20, log="x", cex=.5, mgp=c(2,.5,0),
      ylim=c(0,3.2), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=1.1)
axis (2, seq(0,3), cex.axis=1.1)
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         cty.mns[j] + c(-1,1)*cty.sds[j], lwd=.5)
}
abline(h=ybarbar)
title("ML estimates",cex.main=.9, line=1)



mu0 = mean(param$mu0[1000:nIter])
cty.mns = colMeans(param$mu[1000:nIter, ])
cty.sd = apply(param$mu[1000:nIter, ], 2, sd)

par (mar=c(5,5,4,2)+.1)
plot (sample.size.jittered, cty.mns, cex.lab=1.1, cex.axis=1.1,
  xlab="sample size in county j", ylab=expression (paste
  ("posterior expectation, ", mu[j])),
  pch=20, log="x", cex=.5, mgp=c(2,.5,0), ylim=c(0,3.2), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=1.1)
axis (2, seq(0,3), cex.axis=1.1)
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
    as.vector(cty.mns[j]) + c(-1,1)*cty.sd[j], lwd=.5, col="gray10")
}
abline(h=mu0)
title("Hierarchical Bayesian model",cex.main=.9, line=1)

