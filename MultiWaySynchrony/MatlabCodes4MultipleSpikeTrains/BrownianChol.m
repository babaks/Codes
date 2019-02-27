%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate the cholesky decomposition of the covariance matrix of Brownian Motion
% Input: n is number of time bins
% Output: Cholesky decomposition of the covariance matrix of Brownian Motion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [brownian_chol] = BrownianChol( n )	

brownian_chol = diag(ones(n,1));
for i = 1:(n-1) 
		brownian_chol(i,(i+1):n) = ones( n-i, 1 );
end    
end
