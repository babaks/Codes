% This program performs rejection sampling from Beta(3, 10) using a
% Uniform(0, 1) and M=4.

function rejectionSampling

% Number of samples
n = 3000;

% Ploting the actual distribution for comparision and to see that
% 4*Uniform(0, 1) has the blanket property.
figure;
x0 = [0:0.01:1];
y1 = betapdf(x0, 3, 10);
plot(x0, y1);
y2 = 4*ones(size(x0));
hold on
plot(x0, y2, 'r')
set(gca, 'yLim', [0, 5]);
legend('Beta(3, 10) distribution', 'Uniform(0, 1) multiplied by M = 4')
drawnow
pause(2);
hold off
samp = [];

for i = 1:n
    
    % We first sample from Uniform(0, 1)
    x(i) = unifrnd(0, 1);
    
    % Calculate l(x)/u(x), since the density of uniform is the constant 1,
    % u(x) = Mg(x) = 4;
    p = betapdf(x(i), 3, 10);
    
    % We sample again from Uniform(0, 1) but this time to decide whether we
    % should accept or reject the sample.
    u = unifrnd(0, 4);
    if u<=p
        
        % We accept the sample if u<= p
        samp = [samp; x(i)];
        
        % This part plots the sample and compare it with the actual Beta(3,
        % 10) distribution.
        plot(x0, y1);
        set(gca, 'xLim', [0, 1], 'yLim', [0, 5]);        
        hold on
        plot(x0, y2, 'r')
        [f, c] = hist(samp, [0:.02:1]);
        bar(c, f/(.02*length(samp)), 'FaceColor', 'none');
        legend('Beta(3, 10) distribution', 'Uniform(0, 1) multiplied by M = 4')

        drawnow;
        hold off
    end
end
    