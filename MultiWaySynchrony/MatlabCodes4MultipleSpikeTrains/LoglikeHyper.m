%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the loglikelihood of GP w.r.t. hyperparameters
% Input: yy is the latent variables (n*1 vector), n is number of time bins, logtheta is the logarithm of the hyperparameters and inv_BrownianMat is the invers of the covariance matrix of Brownian Motion
% Output: logpost is the logarithm of the posterior w.r.t. hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function logpost = LoglikeHyper( yy, logtheta, n, inv_BrownianMat )
  
    logprior =  - 0.5*logtheta^2/3^2; % Log-normal priors for hyperparameters
    
    invC = exp(logtheta) * inv_BrownianMat; % Cholesky decomposition of covariance matrix
    logdetC = -n*logtheta + log( DetBrownian(1:n,n)); % Logarithm of the determinant of the covariance matrix
	loglike = -0.5 * logdetC - 0.5 * yy * invC * yy'; % Log-likelihood
	
    logpost = logprior + loglike; 	% log-posterior of hyperparameters

end


    
    
    