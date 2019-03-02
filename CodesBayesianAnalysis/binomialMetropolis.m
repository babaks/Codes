% This program uses the Metropolis algorithm for the binomial model with Beta prior
% Note that this is a conjugate situation so we can compare our results
% with the result based on the closed form.

function theta = binomialtMetropolis
global yObs n alpha beta;

% Number of iterations and number of samples we discard (burn in).
nIterations = 500;
burnIn = 100;

% These are the parameters of the prior Beta(1, 1).
alpha = 1; beta = 1;

% This the observed data
yObs = 39;
n = 100;

% This is our initial point
theta(1) = 0.5;

% This is to create a dynamic plot for MCMC sampling
h = plot(theta);
set(h,'EraseMode','xor','MarkerSize',12)
xlabel(gca, 'Iteration', 'FontSize', 18); ylabel(gca, '\theta', 'FontSize', 18);
axis square
grid off

for i = 2:nIterations

    % Proposing a new point from Uniform(theta-0.1, theta+0.1) distribution
    thetaPrime = unifrnd(theta(i-1) - 0.1, theta(i-1) + 0.1);
    
    % Evaluating the posterior probability of the new point
    postThetaPrime = getPosterior(thetaPrime);
    
    % Evaluating the posterior probability of the current point
    postThetaCurrent = getPosterior(theta(i-1));
    
    % Evaluating the acceptance probability
    acceptanceRate = min(1, postThetaPrime/postThetaCurrent);
    
    % Deciding whether to accept the new point or stay where we are
    u = rand;
    if u<= acceptanceRate
        theta(i) = thetaPrime;
    else
        theta(i) = theta(i-1);
    end
    
    % Ploting the samples
    set(h, 'XData', [1:i], 'YData', theta)
    axis([1 i 0 1])
    drawnow

end


% Comparing the result with what we obtain using the closed form of the posterior distribution 
% based on the conjugacy.
figure
theta0 = [0:0.01:1];
thetaConjugate = betapdf(theta0, 40, 62);
plot(theta0, thetaConjugate);
xlabel(gca, '\theta', 'FontSize', 18); ylabel(gca, 'Density', 'FontSize', 18);
hold on


% [f, c] = hist(theta(burnIn:end), [0:.02:1]);
% bar(c, f/(.02*length(theta(burnIn:end))), 'FaceColor', 'none');

% You can use this code to plot a smooth density
[f, xi] = ksdensity(theta(burnIn:end));
plot(xi, f, 'r')

legend('Conjugate posterior distribution', 'MCMC Posterior distribution')
hold off



% The following function evaluates the psoterior probability for a given
% theta
function posterior = getPosterior(theta)
global yObs n alpha beta;
    
% We use the following to get the posterior distribution. 
% FOR LARGER SAMPLES, IT IS BETTER TO USE THE LOG TRANSFORMATION OF THESE
% DISTRIBUTION AND MODIFY THE ACCEPTANCE PROBABILITY ACCORDINGLY.
prior = theta^(alpha-1)*(1-theta)^(beta-1);
likelihood = theta^(yObs) * (1-theta)^(n-yObs);
posterior = prior*likelihood;

% Alternatively, we could simplify the form of the posterior distribution
% posterior = (theta^(39))*(1-theta)^(61);