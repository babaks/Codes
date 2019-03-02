% This program uses the Metropolis algorithm for the multivariate normal 
% model with unknown mean and known covariance matrix. We use a
% multivariate normal prior for the mean.
% Note that this is a conjugate situation so we can compare our results
% with the result based on the closed form.
% The difference between this program and normPostMetroplis1.m is that this
% oen updates parameters one at a time. 

function theta = normMetropolis2
global yObs n sigma mu0 lambda0;

% Number of iterations and number of samples we discard (burn in).
nIterations = 2000;
% You have to decided on burnin value after you look at the posterior samples.
burnIn = 100;

% These are the parameters of the prior
mu0 = [0, 0];
lambda0 = 10*eye(2);

% This the observed data
yObs = [-1.2, 2.3; -0.5, 0.7; -2.1, -1];
[n, d] = size(yObs);

sigma = cov(yObs);

% This is our initial point
theta(1, :) = [0, 0];

% This is to create a dynamic plot for MCMC sampling

mu_n = (mu0*inv(lambda0) + n*mean(yObs)*inv(sigma))*inv(inv(lambda0)+n*inv(sigma));
lambda_n = inv(inv(lambda0)+n*inv(sigma));
figure
[theta01, theta02] = meshgrid(-4:.2:4,-4:.2:4);
for i = 1:size(theta01, 1) 
for j = 1:size(theta01, 2)
    thetaConjugate(i, j) = mvnpdf([theta01(i, j), theta02(i, j)], mu_n, lambda_n);
end
end
contour(theta01, theta02, thetaConjugate, 'LineWidth', 2)
title('The contour plot shows the conjugate posterior distribution');

xlabel(gca, '\mu_1', 'FontSize', 18); ylabel(gca, '\mu_2', 'FontSize', 18);
hold on
h = plot(theta(1, 1), theta(1, 2), '*', 'MarkerSize', 10);
% set(h,'EraseMode','xor','MarkerSize',12)
xlabel(gca, '\mu_1', 'FontSize', 18); ylabel(gca, '\mu_2', 'FontSize', 18);
axis([-3 0.5 -2 3])


for i = 2:nIterations

    for j = 1:d
    
	thetaPrime = theta(i-1, :);
    % Proposing a new point from N(theta, 0.5) distribution
    thetaPrime(1, j) = normrnd(theta(i-1, j), sqrt(0.5));
    
    % Evaluating the posterior probability of the new point
    postThetaPrime = getPosterior(thetaPrime);
    
    % Evaluating the posterior probability of the current point
    postThetaCurrent = getPosterior(theta(i-1, :));
    
    % Evaluating the acceptance probability
    acceptanceRate = min(1, postThetaPrime/postThetaCurrent);
    
    % Deciding whether to accept the new point or stay where we are
    u = rand;
    if u<= acceptanceRate
        theta(i, j) = thetaPrime(1, j);
    else
        theta(i, j) = theta(i-1, j);
    end


    if i<=30
        if j==1
            plot([theta(i-1, 1), theta(i, 1)], [theta(i-1, 2), theta(i-1, 2)], '--');
        else
            plot([theta(i, 1), theta(i, 1)], [theta(i-1, 2), theta(i, 2)], '--');
            plot(theta(i, 1), theta(i, 2), '*', 'MarkerSize', 10)
        end
        axis square
        drawnow
        pause(7/i)
    end

    
    end

    % Ploting the first 20 samples


end


% Comparing the result to what we obtain using the closed form of the posterior distribution 
% based on the conjugacy.
figure
contour(theta01, theta02, thetaConjugate, 'LineWidth', 2)
axis([-3 0.5 -2 3])
title('The contour plot shows the conjugate posterior distribution');

xlabel(gca, '\mu_1', 'FontSize', 18); ylabel(gca, '\mu_2', 'FontSize', 18);

hold on

plot(theta(burnIn:end, 1), theta(burnIn:end, 2), '.', 'markerSize', 5);
hold off
pause(3)
figure
plot(theta(:, 1), 'b');
hold on
plot(theta(:, 2), 'r');
legend('\mu_1', '\mu_2')
xlabel('Iteration')
title('Trace plot of \mu_1 and \mu_2')
xlabel(gca, 'Iteration', 'FontSize', 18); ylabel(gca, '\mu', 'FontSize', 18);




% The following function evaluates the psoterior probability for a given
% theta
function posterior = getPosterior(theta)
global yObs n sigma mu0 lambda0;
    
% We use the following to get the posterior distribution. 
% FOR LARGER SAMPLES, IT IS BETTER TO USE THE LOG TRANSFORMATION OF THESE
% DISTRIBUTION AND MODIFY THE ACCEPTANCE PROBABILITY ACCORDINGLY.

prior = mvnpdf(theta, mu0, lambda0);
likelihood = prod(mvnpdf(yObs, theta, sigma));
posterior = prior*likelihood;
