% This is a continuous space Markov chain with N(x/2, sqrt(3/4)) transition
% distribution. Note that MATLAB considers the second parameter as standard
% deviation not variance.
function continuousMC1

% Number of samples from the Markov chain
n = 1000;

figure;
set(gca,'xLim', [-4 4], 'yLim', [0, 0.6])
hold on;

% We start at point 0.5
X = 0.5;
for i = 1:n
    % when we are at point X, we sample a new point from N(X/2, sqrt(3/4)).
    X = normrnd(X/2, sqrt(3/4));

    samp(i) = X;

    % This part plots the samples. and shows where the transition
    % distributon is at any given time.
    plot(X, 0, 'rx', 'markersize', 18);
    set(gca,'xLim', [-4 4], 'yLim', [0, 0.6])
    hold on
    x = [X/2-3:.1:X/2+3];
    y = normpdf(x, X/2, sqrt(3/4));
    plot(x, y);
    plot(samp, zeros(size(samp)), 'o')
    drawnow;
   
    title({'Sampling from a Markov chain with N(X/2, 3/4) transition distribution';...
    'The X shows the current stae, and o are all the points we have sampled so far'})    
    hold off
    
    % This is so you can see the first few points slowly.
    pause(4/i);
end
figure
hist(samp(100:end));