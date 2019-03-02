# This program uses logistic regression model to predict the risk of breasst cancer. 
# We are using the Breast Cancer data (bcw.data) from Wisconsin available at UCI repository.
# You need the bcw.data file in the same directory to run this program. 

# these are the step size, w, and the maximum number of stetps, m, we are 
# willing to take
w = 3; 
m = 100;
 
# Number of MCMC iterations, and number of burn-in samples
nIterations = 5000;
burnIn = 500;
 
# Reading in the Wisconsin breast cancer data 
tumor = as.matrix(read.table("bcw.data", sep = ",", header=FALSE));
n = dim(tumor)[1];
dim2 = dim(tumor)[2];
 
# The last column is the outcome variable: 2 for benign, 4 for malignant.
y = rep(0, n); 
y[tumor[, dim2]==4] = 1;
 
# The fist column is the ID number and coulumns 2-10 are the predictors
x = tumor[, 2:(dim2-1)];

p = dim(x)[2];
 
# Standardizing the predictors so they have mean zero and standard
# deviation 1.
x = scale(x);
 
# Randomly choosing 2/3 of data as training and the rest as the test set.
# Note that the result might slightly change depending how you split the
# data. A better approach would be cross-calidation.
ind = sample(seq(1, n), floor(2*n/3))

xTrain = x[ind, ];
yTrain = y[ind];
 
xTest = x[setdiff(seq(1, n), ind), ];
yTest = y[setdiff(seq(1, n), ind)];
 
# Adding a column of 1's for the intercept.
n = dim(xTrain)[1];
xTrain = cbind(rep(1, n), xTrain);
 
# ni is the total number of cases for observation i. In this case, each i
# correspond to only one subject.
ni = rep(1, n);
 
nTest = dim(xTest)[1];
xTest = cbind(rep(1, nTest), xTest);
 
# The parameters of priors for beta's: beta_j ~ N(mu_0j, [tau_0j]^2). Note that all the means are set to
# zzero, and we are using a very broad variance. It is better to regard tau_0 as a hyperparameter for 
# coefficients (excluding the intercept), and give it an inv-chi^2
# hyperprior.
mu0 = rep(0, p+1);
tau02 = 100*rep(1, p+1);
 
# This is the matrix that holds the MCMC samples.
beta = matrix(rep(0, nIterations*(p+1)), nrow = (p+1), ncol=nIterations);

logisticTumor <- function() {
	 
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
postPredProb = matrix(rep(0, nTest*length(burnIn:nIterations)), nrow = nTest, ncol=length(burnIn:nIterations));
logProb = matrix(rep(0, nTest*length(burnIn:nIterations)), nrow = nTest, ncol=length(burnIn:nIterations));
for(i in burnIn:nIterations){
    count = count+1;
    # To get the p(y_tilde | x_tilde), we use the MCMC samples to calculate
    # mu = exp(eta)/(1+exp(eta)), where eta=x_tilde * beta_(s) and beta_(s)
    # is the vectore of beta's at iteration "s".
    eta = xTest%*%beta[, i];    
    postPredProb[, count] = exp(eta)/(1+exp(eta));
    # We also calculate the probability of true value of y_tilde to obtain
    # log-posterior probability at each iteration.
    
    logProb[, count] = log(dbinom(yTest, rep(1, nTest), postPredProb[, count] ));
}
 
# We average over posterior samples to obtain the expected value of posterior predictive
# probabilities
postExp = apply(postPredProb, 1, mean);
 
# Since we are not using any loss function, we predict y_tilde=1 if it's
# posterior probability is bigger than 0.5.
yPred = rep(0, nTest);
yPred[postExp>0.5] = 1;
 
# This give the cross-tabulation of predictive vs. actual values for the
# test cases.
print(table(yTest, yPred))
 
# This calculate the accuracy of our predictions
Accuracy = mean(yTest == yPred)
 
# This calculate the average of log-probabilites over test cases and over
# MCMC samples.
sampLogProb = mean(logProb);
meanLogProb = mean(sampLogProb[is.finite(sampLogProb)])

return(list(beta = beta, Accuracy = Accuracy, postPredProb = postPredProb, meanLogProb = meanLogProb ))
}
 
# This is the slice sampler for regression parameters.
doSliceSampling <- function(currentParam, j, mu0, tau02, allParam){
#  As described by Neal (2003), "In practice, it is often safer to compute 
# g(x) = log(f(x)) rather than f(x)  itself, in order to avoid possible 
# problems with floating-point underflow. One can then use the auxiliary 
# variable z = log(y) = g(x)-e, where e is exponentially distributed with
# mean one, and define the slice by S = {x : z < g(x)}. 
z = getLogPost(currentParam, mu0, tau02, allParam) - rexp(1, 1);

# Stepping out to obtain a [L, R] range  
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
res<-logisticTumor()
