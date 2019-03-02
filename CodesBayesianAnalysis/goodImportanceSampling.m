% This is an example of a good importance sampling. Here we want to use t_3
% for sampling from N(0, 1) so we can approximate the expectation of x^2 
% (with respect to the normal distribution)
function goodImportanceSampling

n = 4000;

% Here, we plot the two distributions
x0 = [-10:0.05:10];
y1 = exp(-0.5*(x0.^2));
hold on
plot(x0, y1)
y2 = (1+(1/3)*x0.^2).^(-2);
plot(x0, y2, '--');
set(gca, 'yLim', [0, 1.1]);
legend('N(0, 1) distribution', 't_3 distribution')

drawnow

pause(2);
hold off
samp = [];
for i = 1:n
    
    % Here, we first sample from t distribjution with 3 degrees of freedom.
    x = trnd(3);
    
    % Then, we calculate the unnormalized density of notmal at the new data
    % point.
    fx = exp(-0.5*(x.^2));
    
    % We also calculate the unnormalized density of t_3 at that point.
    gx = (1+(1/3)*x.^2).^(-2);
    
    % We now calculate the weight
    w(i) = fx/gx;  
    samp(i) = x;
end

% We use the Monte Carlo formula to estimate E(x^2) which would be very close
% to its true value 1.
Ex2 = sum((samp.^2).*w/sum(w))
    