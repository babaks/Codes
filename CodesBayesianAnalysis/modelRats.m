% This is the computer program for the rats tumor model we discussed in the
% class. Here, we use a hierarchical model where we use one parameter
% theta_j for the probability of developing tumor in each rats group.
% Moreover, we assume all theta's share a common Beta(alpha, beta)
% distribution. For simplisity, we fix alpha, and regard beta as a hyperparameter. 
% The hyperprior distribution for beta is Gamma(gamAlpha, gamBeta).
% The outputs of this program are posterior samples for each theta, and the
% posterior samples for beta.

function [theta, beta] = modelRats
global alpha gamAlpha gamBeta

nIterations = 5000;

% Reading in the data.
rats = dlmread('ratTumor.txt');
y = rats(:, 1);
n = rats(:, 2);
nRats = length(y);

% These are the maximum likelihood estimates for theta's.
thetaObs = y./n;

% Setting up the parameters of the priors.
gamAlpha = 2;
gamBeta = 1;
alpha = 2;

% These two matrices hold posterior samples
beta = ones(nIterations, 1);
theta = ones(nIterations, nRats);

for i = 2:nIterations
    % Given the current beta, the posterior distribution of each theta_j
    % has a closed form (Beta distribution) since the prior is conjugate. 
    % So we can use the Gibbs sample to update all the 71 theta's. As you 
    % can see, for each of them, we use their own specific y_j and n_j. 
    for j = 1:nRats
        theta(i, j) = betarnd(alpha+y(j), beta(i-1)+n(j)-y(j));
    end
    
    % Now we move one level higher and take theta's as observations and
    % update beta using a Metropolis step (note, we could use conjugacy and
    % update beta using the Gibbs sampler, but I haven't talked about this
    % in class, so I use Metropolis).
    lower = max(0, beta(i-1)-1);
    upper = beta(i-1)+1;
    
    gCurrent = 1/(upper-lower);
    
    betaPrime = unifrnd(lower, upper);

    lower = max(0, betaPrime-1);
    upper = betaPrime+1;
    
    gPrime = 1/(upper-lower);
        
    % This obtains the log-posterior probabilities for the proposed sample
    % and the current sample.
    logPostPrime = getPost(betaPrime, theta(i, :));
    logPostCurrent = getPost(beta(i-1), theta(i, :));
    
    % Note that I am using the log transformation of posterior
    % distribution, so I am modifying the accceptance probability.
    acceptanceProb = min(1, exp(logPostPrime + log(gPrime) - logPostCurrent - log(gCurrent)));

    % Accepting new proposed values with probability acceptanceProb.
    u = rand;
    if u<= acceptanceProb
        beta(i) = betaPrime;
    else
        beta(i) = beta(i-1);
    end

end

% This part is for creating the plot provided in the course notes
for i = 1:size(theta, 2)
a = prctile(theta(1000:end, i), [2.5, 97.5]);
plot([i, i], [a(1), a(2)])
hold on
xlim([0, i+1])
plot([i], [a(1)+a(2)]/2, 'x')
plot([i+.1], thetaObs(i), 'ro')
end

overallMean = mean(alpha./(alpha+beta(1000:end)));
plot([0, 72], [overallMean, overallMean], 'g')
xlabel(gca, 'Rat Groups', 'FontSize', 14); ylabel(gca, '95% posterior interval for \theta_j', 'FontSize', 14);

% This is the part we calculate the log-posterior probability for a given
% beta and theta. The hyperprior for beta is gamma and the p(theta|alpha,
% beta)=Beta(alpha, beta).
function logPost = getPost(beta, theta)
global alpha gamAlpha gamBeta
logPrior = log(gampdf(beta, gamAlpha, gamBeta));
logLike = sum(log(betapdf(theta, alpha, beta)));
logPost = logPrior + logLike;
