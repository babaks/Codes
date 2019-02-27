function [invSigma] = SolveBrownian( t, n )
invSigma = diag( [ 1./(t(2:n)-t(1:(n-1)))+1./(t(1:(n-1))-[0,t(1:(n-2))]), 1/(t(n)-t(n-1)) ] );
invSigma(1:(n-1),2:n) = invSigma(1:(n-1),2:n) - diag(1./(t(2:n)-t(1:(n-1))));
invSigma(2:n,1:(n-1)) = invSigma(2:n,1:(n-1)) - diag(1./(t(2:n)-t(1:(n-1))));
end
