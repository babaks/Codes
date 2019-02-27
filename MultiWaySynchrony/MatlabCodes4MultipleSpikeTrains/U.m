%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate potential and its gradient w.r.t. beta
% Inputs: beta is the current state of copula model parameters, der is an indicator(0-potential,1-derivative), pmf1smdDATA is used to calculate pmf, dpmfDATA is the derivative of pmf
% Outputs: Potential or its gradient w.r.t. beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[output] = U(beta,der,pmf1smdDATA,dpmfDATA) 

pmf = pmf1smdDATA + dpmfDATA*beta; % Calculate pmf

if der==0
    loglik = sum(log(pmf)); % Loglikelihood
    logpri = -(beta'*beta)/2; % Prior for beta
    output = -(loglik+logpri); % Potential
elseif der==1
    dloglik = sum( dpmfDATA./(pmf*ones(1,length(beta))), 1 ); % Derivative of log-likelihood
%repmat(pmf,1,length(beta))
    dlogpri = -beta';
    output = -(dloglik+dlogpri); % Gradient of potential w.r.t. beta
else
    error('wrong choice of der!');
end

end

