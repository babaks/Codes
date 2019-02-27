%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate the covariance matrix of Brownian Motion
% Input: n is number of time bins
% Output: covariance matrix of Brownian Motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [brownian_mat] = BrownianMat( n )	

brownian_mat = diag(1:n);
for i = 1:(n-1) 
		brownian_mat(i,(i+1):n) = ones( n-i, 1 ) * i;
end
for i = 1:(n-1)
    for j = (i+1):n
			brownian_mat(j,i) = brownian_mat(i,j);
    end
end
    
end
