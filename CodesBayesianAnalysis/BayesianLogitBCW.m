% This program uses logistic regression model to predict the risk of breasst cancer. 
% We are using the Breast Cancer data (bcw.data) from Wisconsin available at UCI repository.
% You need the bcw.data file in the same directory to run this program. 

function [beta, yTest, postPredProb] = BayesianLogitBCW
global xTrain yTrain ni w m

% these are the step size, w, and the maximum number of stetps, m, we are 
% willing to take
w = 3; 
m = 100;

% Number of MCMC iterations, and number of burn-in samples
nIterations = 5000;
burnIn = 500;

% Reading in the Wisconsin breast cancer data 
tumor = dlmread('bcw.data');

% The last column is the outcome variable: 2 for benign, 4 for malignant.
y = logical(tumor(:, end)==4);

% The fist column is the ID number and coulumns 2-10 are the predictors
x = tumor(:, 2:end-1);

[n, p] = size(x);

% Standardizing the predictors so they have mean zero and standard
% deviation 1.
x = x - repmat(mean(x), n, 1);
x = x./repmat(std(x), n, 1);

% Randomly choosing 2/3 of data as training and the rest as the test set.
% Note that the result might slightly change depending how you split the
% data. A better approach would be cross-calidation.
n = size(x, 1);
ind = randsample(1:n, floor(2*n/3));

xTrain = x(ind, :);
yTrain = y(ind, :);

xTest = x(setdiff(1:n, ind), :);
yTest = y(setdiff(1:n, ind), :);

% Adding a column of 1's for the intercept.
n = size(xTrain, 1);
xTrain = [ones(n, 1), xTrain];

% ni is the total number of cases for observation i. In this case, each i
% correspond to only one subject.
ni = ones(n, 1);

nTest = size(xTest, 1);
xTest = [ones(nTest, 1), xTest];

% The parameters of priors for beta's: beta_j ~ N(mu_0j, [tau_0j]^2). Note that all the means are set to
% zzero, and we are using a very broad variance. It is better to regard tau_0 as a hyperparameter for 
% coefficients (excluding the intercept), and give it an inv-chi^2
% hyperprior.
mu0 = zeros(p+1, 1);
tau02 = 100*ones(p+1, 1);

% This is the matrix that holds the MCMC samples.
beta = zeros(p+1, nIterations);

for iter = 2:nIterations
   iter
   % allBeta is the current vector of beta's 
   allBeta = beta(:, iter-1); 
   
   % We cycle through beta's one at a time and update them keeping the
   % other beta's fixed.
   for j = 1:(p+1)
      
       % This is the current beta_j we want to update
       currentBeta = beta(j, iter-1);
       
       % We pass the required information to the slice sampling function,
       % and obtain a new value beta_j at iteration "iter"
       beta(j, iter) = doSliceSampling(currentBeta, j, mu0(j, 1), tau02(j, 1), allBeta);
       
       % We also update allBeta since we have to use the most current
       % values when sampling the nest beta's.
       allBeta(j, 1) = beta(j, iter);
       
   end
end


% We discard some initial samples, and use the rest for predicting the
% outcome for the test cases.

% This part makes the program run slow, you can comment it out if you
% don't need it (i.e., if your objective is not prediction, rather inference 
% bases on posterior distribution of parameters). 
count = 0;
for i = burnIn:nIterations
    count = count+1;
    % To get the p(y_tilde | x_tilde), we use the MCMC samples to calculate
    % mu = exp(eta)/(1+exp(eta)), where eta=x_tilde * beta_(s) and beta_(s)
    % is the vectore of beta's at iteration "s".
    eta = xTest*beta(:, i);    
    postPredProb(:, count) = exp(eta)./(1+exp(eta));
    % We also calculate the probability of true value of y_tilde to obtain
    % log-posterior probability at each iteration.
    
    logProb(:, count) = log(binopdf(yTest, ones(nTest, 1), postPredProb(:, count) ));
end

% We average over posterior samples to obtain the expected value of posterior predictive
% probabilities
postExp = mean(postPredProb, 2);

% Since we are not using any loss function, we predict y_tilde=1 if it's
% posterior probability is bigger than 0.5.
yPred = logical(postExp>0.5);

% This give the cross-tabulation of predictive vs. actual values for the
% test cases.
crosstab(yTest, yPred)

% This calculate the accuracy of our predictions
Accuracy = mean(yTest == yPred)

% This calculate the average of log-probabilites over test cases and over
% MCMC samples.
sampLogProb = mean(logProb);
meanLogProb = mean(sampLogProb(~isinf(sampLogProb)))


% This is the slice sampler for regression parameters.
function newParam = doSliceSampling(currentParam, j, mu0, tau02, allParam)
global w m;

%  As described by Neal (2003), "In practice, it is often safer to compute 
% g(x) = log(f (x)) rather than f(x)  itself, in order to avoid possible 
% problems with floating-point underflow. One can then use the auxiliary 
% variable z = log(y) = g(x) -e, where e is exponentially distributed with
% mean one, and define the slice by S = {x : z < g(x)}."
 
z = getLogPost(currentParam, mu0, tau02, allParam) - exprnd(1);

% Stepping out to obtain the [L, R] range
u = rand;
L = currentParam - w*u;
R = L + w;
v = rand;
J = floor(m*v);
K = (m-1) - J;

allParam(j) = L;
while J>0 && z < (getLogPost(L,  mu0, tau02, allParam))
    L = L - w;
    J = J - 1;
    allParam(j) = L;
end

allParam(j) = R;
while K>0 && z < (getLogPost(R,  mu0, tau02, allParam))
    R = R+w;
    K = K-1;
    allParam(j) = R;
end


% Shrinkage to obtain a sample
u = rand;
newParam = L + u*(R-L);
allParam(j) = newParam;
while z > (getLogPost(newParam, mu0, tau02, allParam))
    if newParam < currentParam
        L = newParam;
    else
        R = newParam;
    end
    
    u = rand;
    newParam = L + u*(R-L);
    allParam(j) = newParam;
end



% This function calculates the log-posterior probability based on the logit
% model with a normal prior for beta_j.
function logPosterior = getLogPost(currentBeta, mu0, tau02, allBeta)
global xTrain yTrain ni;

eta = xTrain*allBeta;

logLike = sum(yTrain.*eta - ni.*log(1+exp(eta)));

logPrior =   -( (currentBeta-mu0).^2 ) / (2*tau02);

logPosterior = logLike + logPrior;

