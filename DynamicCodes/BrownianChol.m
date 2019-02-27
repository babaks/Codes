function [brownian_chol] = BrownianChol( n )	

brownian_chol = diag(ones(n,1));
for i = 1:(n-1) 
		brownian_chol(i,(i+1):n) = ones( n-i, 1 );
end    
end
