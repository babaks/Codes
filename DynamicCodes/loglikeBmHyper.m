function logpost = loglikeBmHyper( yy, logtheta, n, invBrownianMat )
    %function to calculate the loglikelihood of GP w.r.t. hyperparameters
    % yy is the outputs (n*1 vector)
    % n is number of outputs
    % DiffMat is the matrix to calculate the covariance matrix of the
    % Gaussian Process
    
    % Log-normal priors for hyperparameters
    logprior =  - 0.5*logtheta^2/3^2;
    
    % Covariance Matrix
    invC = exp(logtheta) * invBrownianMat; % Cholesky decomposition
    logdetC = -n*logtheta + log( DetBrownian(1:n,n)); % Logarithm of the determinant of the covariance matrix
	loglike = -0.5 * logdetC - 0.5 * yy * invC * yy';

	% log-posterior of hyperparameters
	logpost = logprior + loglike;
end
