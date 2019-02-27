function logpost = loglikeGpHyper( yy, logeta, logrho, logalpha, logJ2, n, DiffMat )
    %function to calculate the loglikelihood of GP w.r.t. hyperparameters
    % yy is the outputs (n*1 vector)
    % n is number of outputs
    % DiffMat is the matrix to calculate the covariance matrix of the
    % Gaussian Process
    
    % Log-normal priors for hyperparameters
    logprior =  -0.5*logrho^2/3^2 - 0.5*logeta^2/3^2 - 0.5*logalpha^2/3^2 - 0.5*logJ2^2/3^2;
    
    % Covariance Matrix
    CC = exp(logJ2) * diag(ones(1,n)) + exp(logeta) * exp( -DiffMat * exp(logrho) ) + exp(logalpha)*ones(n);
    LL = chol( CC ); % Cholesky decomposition
    logdetC = 2 * sum( log(diag(LL)) ); % Logarithm of the determinant of the covariance matrix
    invLL = inv( LL ); % Inverse of LL
	invLLy = transpose(invLL) * yy;
	loglike = -0.5 * logdetC - 0.5 * transpose(invLLy) * invLLy;

	% log-posterior of hyperparameters
	logpost = logprior + loglike;
end


    
    
    