% This program uses the Gibbs sampler for simulating from the posterior
% distribution of a normal model with unknown mean and variance. The mean
% and percision (i.e., 1/variance) are given normal and Gamma priors. 

function [mu, sigma2] = normPostGibbs
global yObs n mu0 sigma02 alpha beta;

nIterations = 200;

mu0 = 65;
sigma02 = 100^2;

alpha = 1;
beta = 1;

yObs = [72, 75, 70]';
[n, d] = size(yObs);

mu(1) = 50;
sigma2(1) = 1;

% This is to create a dynamic plot for MCMC sampling
h1 = subplot(2, 1, 1);
plot(mu);

xlabel(h1, '\iteration', 'FontSize', 18); ylabel(h1, '\mu', 'FontSize', 18);
axis tight
grid off

h2 = subplot(2, 1, 1);
plot(1, sqrt(sigma2));

xlabel(h2, '\iteration', 'FontSize', 18); ylabel(h2, '\sigma', 'FontSize', 18);
axis tight
grid off

for i = 2:nIterations

    % The conditional distribution of mu given sigma is normal(mu_n,
    % sigma_n) when we use a normal(mu0, sigma02) prior for mu.
    mu_n = (mu0/sigma02 + sum(yObs)/sigma2(i-1)) / (1/sigma02 + n/sigma2(i-1));
    sigma_n = sqrt(1/(1/sigma02 + n/sigma2(i-1)));
    mu(i) = normrnd(mu_n, sigma_n);
    
    % The conditional distribution of precision (i.e., 1/sigma2) is Gamma 
    % when we use a Gamma(alpha, beta) prior. Note that if you are using R,
    % you have to inverse beta, and also inverse the second parameter of the
    % posterior distribution.
    precision = gamrnd(alpha + n/2, 1/(1/beta + sum( (yObs - mu(i)).^2 )/2));
    sigma2(i) = 1/precision;
    
    % Ploting the samples
    h1 = subplot(2, 1, 1);    
    plot(mu)
    xlabel(h1, 'iteration', 'FontSize', 18); ylabel(h1, '\mu', 'FontSize', 18);
   
    h2 = subplot(2, 1, 2);    
    plot(sqrt(sigma2))

    xlabel(h2, 'iteration', 'FontSize', 18); ylabel(h2, '\sigma', 'FontSize', 18);

    drawnow
        
end



