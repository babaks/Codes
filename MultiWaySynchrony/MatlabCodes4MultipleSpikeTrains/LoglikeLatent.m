%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate log-likelihood w.r.t. latent variables
% Inputs: u is the latent variables and y is the data
% Outputs: output is the log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = LoglikeLatent(u,y) 

[N,n] = size(y); % N is number of trials and n is number of time bins
p = exp(u) ./ (1+exp(u)); % Calculate firing probabilities
loglik = sum(sum(y .* log(ones(N,1)*p) + (1-y) .* log(1-ones(N,1)*p))); % Calculate log-likelihood
output = loglik;

end
