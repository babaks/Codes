% This is an example where the importance sampling will not work properly.
% Here, we are using N(0, 1) to sample from t_3.

function badImportanceSampling

n = 1000;

% This part plots N(0, 1) and t_3
hold off
x0 = [-10:0.05:10];
y1 = (1+(1/3)*x0.^2).^(-2);
plot(x0, y1);
y2 = exp(-0.5*(x0.^2));
hold on
plot(x0, y2, '--')
set(gca, 'yLim', [0, 1.5]);
drawnow
legend('t_3 distribution', 'N(0, 1) distribution')
pause(2);
hold off
samp = [];
for i = 1:n
    
    % We sample from N(0, 1)
    x = normrnd(0, 1);
    
    % Calculate the unnormalized density of t_3
    fx = (1+(1/3)*x.^2).^(-2);
    
    % Calculate the unnormalized density of N(0, 1)
    gx = exp(-0.5*(x.^2));
    
    % Calculate the weight
    w(i) = fx/gx;  
    samp(i) = x;
end

% Use the Monte Carlo formula to estimate the expectation of x^2 with
% respect to t_3 distribution. The answer would be systematically lower
% than the actual value 3.
Ex2 = sum((samp.^2).*w/sum(w))
    