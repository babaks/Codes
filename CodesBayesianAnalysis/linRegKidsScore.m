% This is a Bayesian linear regression model with informative priors
% The outputs are MCMC samples for the regression parameters "beta", the
% variance of y|x.beta (sigma2). Note that we are using a simple Gibbs sampler here.
function [beta, sigma2] = linRegKidsScore
global n d nu0 sigma02 mu0 tau02 nIterations

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

% Here, we are standardizing predictors by subtracting mean and deviding by
% standard deviation.
x = x - repmat(mean(x), n, 1);
x = x ./ repmat(std(x), n, 1);


% d is the number of regression parameters including the intercept.
d = p+1;
x = [ones(n, 1), x];

% The parameteres of Inv-chi2 prior for sigma2
nu0 = 1;
sigma02 = 0.5;

% The parameters of normal prior N(mu0, tau02) for beta's.
mu0 = zeros(1, d);
tau02 = 1000*ones(1, d);

% If n is bigger than d, the posterior sampling is done faster and it 
% is more straightforward.
if n > d
     [beta, sigma2] = regLargeN(y, x);
else
     [beta, sigma2] = regSmallN(y, x);
end



function [beta, sigma2] = regLargeN(y, x)
global n d nu0 sigma02 mu0 tau02 nIterations

% The following two matrices hold the posterior samples for sigma2 and
% beta.
sigma2 = ones(nIterations, 1);
beta = ones(d, nIterations);

Lambda0 = diag(tau02);
invLambda0 = diag(1./tau02);


% It is always better to use the Cholosky decomposition for getting
% the iverse of matrices and sampling from multivariate normal.    
% Here, we first get, L1, the Cholesky decomposition of (x'x). To get, the
% inv(x'x), we use inv(L1)*(inv(L1))'.

X = x'*x;
L1 = chol(X);
invL1 = inv(L1);
invX = invL1*invL1';

% This give the mean of the posterior distribution.
betaHat = invX*x'*y; 

% Since later on we need the sample from N(betaHat, invX*sigma2), we should
% find the Cholesky decomposition of invX*sigma2. Instead of doing this in
% a loop, we can get the Cholesky of invX (we call it L2) and then at each
% iteration, we multiply this by sqrt(sigma2).
L2 = chol(invX);

for i = 2:nIterations
    
    % The posterior distribution of beta in this case is N(mu_n,
    % Lambda_n). Note that this is the same as what we had for multivariate
    % normal model with conjugate priors. However, there is no "n" here
    % since it is already summed over all data points.
   
    L3 = L2*sqrt(sigma2(i-1));
    invL3 = inv(L3);
    invSigma_n = invL3*invL3';
    
    invLambda_n = invLambda0 + invSigma_n;
    L4 = chol(invLambda_n);
    invL4 = inv(L4);
    Lambda_n = invL4*invL4';
    
    mu_n = Lambda_n*(invLambda0*mu0' + invSigma_n*betaHat);
            
    
    % We can now sample from the posterior distribution beta given sigma^2. As
    % discussed in the class, we use the Cholesky decomposition for this.
    
    L5 = chol(Lambda_n);
    u = normrnd(0, 1, 1, d); 
    beta(:, i) = (u*L5)' + mu_n;

    eta = x*beta(:, i);
    eps = y - eta;
    nu = sum(eps.^2)/(n);
    nu_n = nu0 + n;
    sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
    
    % I sample a new value for sigma2 from Inv-chi2(nu_n, sigma02_n), to do
    % this, I sample from z ~ chi2(nu_n) and then my sample from Inv-chi2
    % would be nu_n*sigma02_n/z.
    z = chi2rnd(nu_n);
    sigma2(i) = nu_n*sigma02_n/z;
    
    
end





function [beta, sigma2] = regSmallN(y, x)
global n d nu0 sigma02 mu0 tau02 nIterations

% The following two matrices hold the posterior samples for sigma2 and
% beta.
sigma2 = ones(nIterations, 1);
beta = ones(d, nIterations);

Lambda0 = diag(tau02);

% These are y*, x* 
yStar = [y; mu0'];
xStar = [x; eye(d)];

sigmaStar = blkdiag(eye(n), Lambda0);

for i = 2:nIterations
    
    %  A contains the diagonal elements of SigmaStar^(-1/2)
    A = 1./sqrt([sigma2(i-1)*ones(1, n), tau02]');
    % newX is SigmaStar^(-1/2)*xStar
    newX = repmat(A, 1, d).*xStar;
    %  newY is SigmaStar^(-1/2)*yStar
    newY = A.*yStar;
    
    % Now we can get mu_n and Lambda_n as described in the course note based
    % on newX and newY
    
    L2 = chol(newX'*newX);
    invL2 = inv(L2);
    Lambda_n = invL2*invL2';
    mu_n = Lambda_n*newX'*newY;
    
    % We can now sample from the posterior distribution of beta given sigma^2. As
    % discussed in the class, we use the Cholesky decomposition for this.
    L3 = chol(Lambda_n);
    u = normrnd(0, 1, 1, d);
    beta(:, i) = (u*L3)' + mu_n;
    
    % Now, given the current beta's, we sample from the posterior
    % distribution of sigma2, which is Inv-chi2(nu_n, sigma02_n).
    eta = x*beta(:, i);
    eps = y - eta;
    nu = sum(eps.^2)/(n);
    nu_n = nu0 + n;
    sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
    
    % I sample a new value for sigma2 from Inv-chi2(nu_n, sigma02_n), to do
    % this, I sample from z ~ chi2(nu_n) and then my sample from Inv-chi2
    % would be nu_n*sigma02_n/z.
    z = chi2rnd(nu_n);
    sigma2(i) = nu_n*sigma02_n/z;
    
    
end


