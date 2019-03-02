% This program uses the Metropolis-Hastings (MH) algorithm for the Poisson 
% model with Gamma prior Note that this is a conjugate situation so we can 
% compare our results with the result based on the closed form.

function theta = poissonMH
global yObs n alpha beta;

% Number of iterations and number of samples we discard (burn in).
nIterations = 1000;
burnIn = 200;

% These are the parameters of the prior
alpha = 3; beta = 1/15;

% This the observed data
yObs = [0, 1];
n = length(yObs);

% This is our initial point
theta(1) = 0.5;

% This is to create a dynamic plot for MCMC sampling
h = plot(theta);
set(h,'EraseMode','xor','MarkerSize',12)
xlabel(gca, 'Iteration', 'FontSize', 18); ylabel(gca, '\theta', 'FontSize', 18);
axis square
grid off

for i = 2:nIterations

    % Proposing a new point from the log-normal distribution distribution
    thetaPrime = lognrnd(log(theta(i-1)), sqrt(0.5));
    
    % Evaluating the posterior probability of the new point
    postThetaPrime = getPosterior(thetaPrime);
    
    % Evaluating the posterior probability of the current point
    postThetaCurrent = getPosterior(theta(i-1));
    
    % Transition probability from the current theta to theta'
    g_thetaCurrent_ThetaPrime = lognpdf(thetaPrime, log(theta(i-1)), sqrt(0.5));
    
    % Transition probability from theta' to the current theta
    g_thetaPrime_thetaCurrent = lognpdf(theta(i-1), log(thetaPrime), sqrt(0.5));
    
    % Evaluating the acceptance probability
    acceptanceRate = min(1, (postThetaPrime*g_thetaPrime_thetaCurrent)/(postThetaCurrent*g_thetaCurrent_ThetaPrime));
    
    % Deciding whether to accept the new point or stay where we are
    u = rand;
    if u<= acceptanceRate
        theta(i) = thetaPrime;
    else
        theta(i) = theta(i-1);
    end
    
    % Ploting the samples
    set(h, 'XData', [1:i], 'YData', theta)
    axis([1 i 0 max(theta)+.2])
    drawnow

end


% Comparing the result with what we obtain using the closed form of the posterior distribution 
% based on the conjugacy.
figure
theta0 = [0:0.01:2];
thetaConjugate = gampdf(theta0, 4, 1/17);
plot(theta0, thetaConjugate);
xlabel(gca, '\theta', 'FontSize', 18); ylabel(gca, 'Density', 'FontSize', 18);
hold on


% [f, c] = hist(theta(burnIn:end), [0:.02:1]);
% bar(c, f/(.02*length(theta(burnIn:end))), 'FaceColor', 'none');

% You can use this code to plot a smooth density
[f, xi] = ksdensity(theta(burnIn:end));
plot(xi(xi>0), f(xi>0), 'r')


legend('Conjugate posterior distribution', 'MCMC Posterior distribution')
hold off


% The following function evaluates the psoterior probability for a given
% theta
function posterior = getPosterior(theta)
global yObs n alpha beta;
    
% We use the following to get the posterior distribution. 
% FOR LARGER SAMPLES, IT IS BETTER TO USE THE LOG TRANSFORMATION OF THESE
% DISTRIBUTION AND MODIFY THE ACCEPTANCE PROBABILITY ACCORDINGLY.
prior = theta^(alpha-1)*exp(-theta/beta);
likelihood = theta^(sum(yObs))*exp(-n*theta);
posterior = prior*likelihood;

% Alternatively, we could simplify the form of the posterior distribution
% posterior = (theta^3)*exp(-17*theta);