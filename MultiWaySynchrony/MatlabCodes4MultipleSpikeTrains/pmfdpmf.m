%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate pmf and derivative of pmf for spike train in 1 bin
% Inputs: y is the data (D*1 vector) and p is the firing probabilities
% Outputs: pmf1smd, which is not pmf but will be used to calculate pmf, and derivative of pmf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = pmfdpmf(y,p)

q = 1-p; % Probabilities of absence of spikes
D = length(y); % Number of neurons
pmf1smd = zeros(1); % Initialization for pmf1smd, which will be used to calculate pmf later
dpmf = zeros(1, D*(D-1)/2 ); % Initialization derivative of pmf

%% Calculate pmf and dpmf
Ones = sum(y); % Number of spikes
M = 2^Ones; 
y_prm = zeros(M,D);
if Ones > 1
	y_prm(:,logical(y)) = combrep(0:1,Ones);
end
if Ones == 1
    y_prm(:,logical(y)) = [0,1];
end
y_diff = ones(M,1) * y - y_prm;
    
mgcdf = (ones(M,1) * q).^(1-y_diff); % marginal cdf, (M,D) matrix
prod2 = prod(mgcdf,2); % Product of marginal cdf
prod3 = zeros(M,D*(D-1)/2);
for i=1:M
	B = (1-mgcdf(i,:))'*(1-mgcdf(i,:));
	prod3(i,:) = B(triu(true(D),1));
end
djtcdf = prod3.*(prod2*ones(1,D*(D-1)/2)); %  Derivative of joint cdf, (M, D(D-1)/2) matrix
sgn = (-1).^sum(y_prm,2);
dpmf = sum((sgn*ones(1,D*(D-1)/2)).*djtcdf,1); % Derivative of pmf, !emphasize sum over column, otherwise get scalar for row vector 
pmf1smd = sum(sgn.*prod2);
output = [pmf1smd,dpmf];
end
