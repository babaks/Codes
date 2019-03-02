% This is a Bayesian linear regression model with noninformative priors
% The outputs are MCMC samples for the regression parameters "beta" and sigma2, the
% variance of (y|x, beta). We are using a simple Gibbs sampler here.
% 
function [beta, sigma2] = linRegKidsScoreNonInfo

% Number of MCMC iterations
nIterations = 10000;

% Reading in the kids test scores dataset
data = importdata('kidiq.txt');
kidiq = data.data;
kts = kidiq(:, 1);
momHs = kidiq(:, 2);
momIq= kidiq(:, 3);

% We choose mother's IQ and whether she is graduated from high school as
% predictors
x = [momHs, momIq];

% The output is kid's test score
y = kts;

[n, p] = size(x);

% Here, we are standardizing the predictors by subtracting the mean and deviding by
% their standard deviation.
x = x - repmat(mean(x), n, 1);
x = x ./ repmat(std(x), n, 1);

% Adding a column of 1's for the intercept
x = [ones(n, 1), x];

% Number of parameters including the intercept.
d = p+1;

% The following two matrices hold the posterior samples for sigma2 and
% beta.
sigma2 = ones(nIterations, 1);
beta = ones(d, nIterations);

% It is always better to use the Cholosky decomposition to invert
% matrices and sampling from multivariate normal.    
% Here, we first get, L1, the Cholesky decomposition of (x'x). To get, the
% inv(x'x), we use inv(L1)*(inv(L1))'.

X = x'*x;
L1 = chol(X);
invL1 = inv(L1);
invX = invL1*invL1';

% This give the mean of the posterior distribution.
betaHat = invX*x'*y;

% Since later on we need the sample from N(betaHat, invX*sigma2), we need
% find the Cholesky decomposition of invX*sigma2. Instead of doing this in
% a loop, we can get the Cholesky of invX (we call it L2) and then at each
% iteration, we multiply this by sqrt(sigma2).
L2 = chol(invX);

for i = 2:nIterations
    
    L3 = L2*sqrt(sigma2(i-1));       
    
    % We can now sample from the posterior distribution of beta given sigma^2. As
    % discussed in the class, we use the Cholesky decomposition for this.
    z = normrnd(0, 1, 1, d);
    beta(:, i) = (z*L3)' + betaHat;
    

    % Now, given the current beta's, we sample from the posterior
    % distribution of sigma2, which is Inv-chi2(n-p-1, s2).
    eta = x*beta(:, i);
    eps = y - eta;
    s2 = sum(eps.^2)/(n-d);

    % I sample a new value for sigma2 from Inv-chi2(n-p-1, s2). To do
    % this, I sample from z ~ chi2(n-p-1) and then my sample from Inv-chi2
    % would be (n-p-1)*s2/z.    
    z = chi2rnd(n-d);
    sigma2(i) = (n-d)*s2/z;
    
    
end
