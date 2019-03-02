# This program uses logistic regression to model the effect of snoring on heat disease.
 
# these are the step size, w, and the maximum number of stetps, m, we are 
# willing to take
w = 3; 
m = 100;
 
# Number of MCMC iterations, and number of burn-in samples
nIterations = 5000;
burnIn = 500;
 
snoreData = as.matrix(rbind(c(0, 24, 1379), c(2, 35, 638),  c(4, 21, 213),  c(5, 30, 254)));

xTrain = snoreData[, 1];
yTrain = snoreData[, 2];
ni = snoreData[, 3];

n = length(xTrain)
p = 1;
 
# Adding a column of 1's for the intercept.
xTrain = cbind(rep(1, n), xTrain);
 
# The parameters of priors for beta's: beta_j ~ N(mu_0j, [tau_0j]^2). Note that all the means are set to
# zzero, and we are using a very broad variance. It is better to regard tau_0 as a hyperparameter for 
# coefficients (excluding the intercept), and give it an inv-chi^2
# hyperprior.
mu0 = rep(0, p+1);
tau02 = 100*rep(1, p+1);
 
# This is the matrix that holds the MCMC samples.
beta = matrix(rep(0, nIterations*(p+1)), nrow = (p+1), ncol=nIterations);

logisticSnore <- function() {
	 
for(iter in 2:nIterations) {
   iter
   # allBeta is the current vector of beta's 
   allBeta = beta[, iter-1]; 
   
   # We cycle through beta's one at a time and update them keeping the
   # other beta's fixed.
   for(j in 1:(p+1)) {
      
       # This is the current beta_j we want to update
       currentBeta = beta[j, iter-1];
       
       # We pass the required information to the slice sampling function,
       # and obtain a new value beta_j at iteration "iter"
       beta[j, iter] = doSliceSampling(currentBeta, j, mu0[j], tau02[j], allBeta);
       
       # We also update allBeta since we have to use the most current
       # values when sampling the nest beta's.
       allBeta[j] = beta[j, iter];
       
   }
}
 

 
# We discard some initial samples, and use the rest for predicting the
# outcome for the test cases.
 
# This part makes the program run slow, you can comment it out if you
# don't need it (i.e., if your objective is not prediction, rather inference 
# bases on posterior distribution of parameters). 
count = 0;
estMu = matrix(rep(0, n*length(burnIn:nIterations)), nrow = n, ncol=length(burnIn:nIterations));
estLogOdds = matrix(rep(0, n*length(burnIn:nIterations)), nrow = n, ncol=length(burnIn:nIterations));
for(i in burnIn:nIterations){
    count = count+1;
    # To get the p(y_tilde | x_tilde), we use the MCMC samples to calculate
    # mu = exp(eta)/(1+exp(eta)), where eta=x_tilde * beta_(s) and beta_(s)
    # is the vectore of beta's at iteration "s".
    eta = xTrain%*%beta[, i];    
 	estMu[, count] = exp(eta)/(1+exp(eta));
    estLogOdds[, count] = eta; 

}
 
postExp = apply(estMu, 1, mean);
postExpLogOdds = apply(estLogOdds, 1, mean);
 
return(list(beta = beta, postExp = postExp, postExpLogOdds = postExpLogOdds ))
}
 
# This is the slice sampler for regression parameters.
doSliceSampling <- function(currentParam, j, mu0, tau02, allParam){
#  As described by Neal (2003), "In practice, it is often safer to compute 
# g(x) = log(f (x)) rather than f(x)  itself, in order to avoid possible 
# problems with floating-point underflow. One can then use the auxiliary 
# variable z = log(y) = g(x) -e, where e is exponentially distributed with
# mean one, and define the slice by S = {x : z < g(x)}." 
z = getLogPost(currentParam, mu0, tau02, allParam) - rexp(1, 1);

# Stepping out to obtain the [L, R] range  
u = runif(1);
L = currentParam - w*u;
R = L + w;
v = runif(1);
J = floor(m*v);
K = (m-1) - J;
 
allParam[j] = L;
while(J>0 && z < (getLogPost(L,  mu0, tau02, allParam)) ){
    L = L - w;
    J = J - 1;
    allParam[j] = L;
}
 
allParam[j] = R;
while(K>0 && z < (getLogPost(R,  mu0, tau02, allParam))){
    R = R+w;
    K = K-1;
    allParam[j] = R;
}

# Shrinkage to obtain a new sample
u = runif(1);
newParam = L + u*(R-L);
allParam[j] = newParam;
while(z > (getLogPost(newParam, mu0, tau02, allParam))){
    if(newParam < currentParam){
        L = newParam;
    }else{
        R = newParam;
    }
    
    u = runif(1);
    newParam = L + u*(R-L);
    allParam[j] = newParam;
}
 
return(newParam) 
} 
 
# This function calculates the log-posterior probability based on the logit
# model with a normal prior for beta_j.
getLogPost <- function(currentBeta, mu0, tau02, allBeta){
 
eta = xTrain%*%allBeta;
 
logLike = sum(yTrain*eta - ni*log(1+exp(eta)));
 
logPrior =   -( (currentBeta-mu0)^2 ) / (2*tau02);
 
logPosterior = logLike + logPrior;

return(logPosterior)

}





# Running the program
res<-logisticSnore()
