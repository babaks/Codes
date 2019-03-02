% This program uses logistic regression to model the effect of snoring on heat disease.

function [beta, postExpMu, postExpLogOdds] = BayesianLogitSnore
global xTrain yTrain ni w m

% these are the step size, w, and the maximum number of stetps, m, we are 
% willing to take
w = 3; 
m = 10;

% Number of MCMC iterations, and number of burn-in samples
nIterations = 1500;
burnIn = 500;
snoreData = [0 24 1379; 2 35 638; 4 21 213; 5 30 254];

xTrain = snoreData(:, 1);
yTrain = snoreData(:, 2);
ni = snoreData(:, 3);
n = size(yTrain, 1);

% Adding a column of 1's for the intercept.
xTrain = [ones(n, 1), xTrain];
[n, p] = size(xTrain);

% The parameters of priors for beta's: beta_j ~ N(mu_0j, [tau_0j]^2). Note that all the means are set to
% zzero, and we are using a very broad variance. It is better to regard tau_0 as a hyperparameter for 
% coefficients (excluding the intercept), and give it an inv-chi^2
% hyperprior.
mu0 = zeros(p, 1);
tau02 = 100*ones(p, 1);

% This is the matrix that holds the MCMC samples.
beta = zeros(p, nIterations);

for iter = 2:nIterations
   iter

   % allBeta is the current vector of beta's 
   allBeta = beta(:, iter-1); 

   % We cycle through beta's one at a time and update them keeping the
   % other beta's fixed.   
   for j = 1:(p)
      
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

% These the heart disease rate and odds
obsMu = snoreData(:, 2)./snoreData(:, 3);
obsLogOdds = log(obsMu./(1- obsMu));

% We discard some initial samples, and use the rest for predicting the
% outcome for the test cases.
count = 0;
for i = burnIn:nIterations
    count = count+1;
    eta = xTrain*beta(:, i);
    estMu(:, count) = exp(eta)./(1+exp(eta));
    estLogOdds(:, count) = eta; 
end

% We average over posterior samples to obtain the expected value of posterior predictive
% probabilities
postExpMu = mean(estMu, 2);
postExpLogOdds = mean(estLogOdds, 2);


% plot(snoreData(:, 1), obsMu, 'o')
% hold on
% plot(snoreData(:, 1), postExpMu, '-- *')
% xlabel('Soring level', 'fontSize', 15)
% ylabel('\mu', 'fontSize', 15)
% legend('Observed proportions', 'posterior predictive probabilities')
% xlim([0, 6])
% 
% figure
% plot(snoreData(:, 1), obsLogOdds, 'o')
% hold on
% plot(snoreData(:, 1), postExpLogOdds, '--*')
% xlabel('Soring level', 'fontSize', 15)
% ylabel('log[\mu/ (1-\mu)]', 'fontSize', 15)
% xlim([0, 6])
% legend('Observed log-odds', 'posterior predictive log-odds')
% 




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



function logPosterior = getLogPost(currentBeta, mu0, tau02, allBeta)
global xTrain yTrain ni;

eta = xTrain*allBeta;

logLike = sum(yTrain.*eta - ni.*log(1+exp(eta)));

logPrior =   -( (currentBeta-mu0).^2 ) / (2*tau02);

logPosterior = logLike + logPrior;
