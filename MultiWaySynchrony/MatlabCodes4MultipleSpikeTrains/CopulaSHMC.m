%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to do mcmc sampling for latent variables, hyperparameters of Brownian Motion, and parameters of the copula model using elliptical slice 
% sampler, slice sampler and spherical HMC, respectively.
% Input: Y is the data (Spike Trains), Niteration and BurnIn are the number of iterations and burn-in of MCMC, w and m are step size and step number of slice sampler.
% Output: Logtheta is the log of the hyperparameter of Brownian Motion, PP is the sampledd marginal firing probabilities, Samp is the sample of copula parameters, 
% logwt is the log weights of the spherical HMC samples, acpt is the acceptance rate of spherical HMC, and TIME is the time to do MCMC sampling.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Logtheta, PP, Samp, logwt, acpt, TIME] = CopulaSHMC( Y, Niteration, BurnIn, w, m )
	
    [N,n,D] = size(Y); % dimension of data; N: number of trials; n: number of time bins; D: number of neurons;
	
    %% Initialization for hyperparameters and latent variables
	Logtheta = zeros( D, Niteration )+3;  % Log of the hyperparameters
    UU = normrnd( -4, 0.1, [D,n] ); % Latent variables
    PP = zeros( D, n, Niteration ); % Marginal firing probabilities
    %% Brownian Motion
    inv_BrownianMat = SolveBrownian(1:n,n); % Inverse of the covariance matrix of Brownian Motion
    brownian_chol = BrownianChol(n); % Cholesky Decomposition of the covariance matrix of Brownian Motion
    
	%% Spherical HMC setting
	TrjLength = 0.05*2*pi/D;   % Length of Trajectory	
	NLeap = 10; % Number of leap frog
	Stepsz = TrjLength/NLeap; % Step size
	%% Storage for spherical HMC
	K = D*(D-1)/2; % Number of parameters of copula model
	Samp = zeros( K, Niteration-BurnIn ); % Variable to store samples of copula model parameters
    acpt = 0; % overall acceptance
	accp = 0; % online acceptance
	logwt = zeros( 1, Niteration-BurnIn ); % Log weight of spherical samples    
	%% Initialization for spherical HMC
	beta = zeros( K, 1 ); % Parameters of copula model
	Beta = zeros( K, Niteration ); % Variable to store the samples of copula model parameters
	theta = times( sign(beta), sqrt(abs(beta)) ); % must have |theta|<=1!
	theta = [theta',sqrt(abs(1-sum(times(theta,theta))))]'; % Map the copula model parameters to the standard sphere

    
    tic;  % Time to start MCMC 
	for i = 1:(Niteration-1)
        
        %% display online acceptance rate per 100 iteration and samples of
        %% copula model parameters and latent variables
        if mod(i,100)==0
            disp( [ 'Iteration ', num2str(i), ' completed!' ] );
            disp( ['Acceptance Rate of Latent Variables: ', num2str(accp/100) ] );
            accp = 0;
        end         
        disp( i );     
        disp( beta(1:5) );
        disp( UU(1,1:5) );

       %% Sample Latent Variables
        % tic;
        parfor j=1:D % Use parfor to sample latent variables in parallel
            ytemp = zeros(N,n);
            ytemp(:,:) = Y(:,:,j); % data of the j^th neuron
			cholC = sqrt(1/exp(Logtheta(j,i))) * brownian_chol; % Cholesky Decoposition of the covariance matrix
            UU(j,:) = SamplerLatent( UU(j,:), ytemp, n, cholC ); % Sample latent variables by calling the elliptical slice sampler.
        end
        PP(:,:,i+1) = ones(D,n) - 1 ./ (1+exp(UU)); % save sampled marginal firing probabilities.   
        % toc
       
        %% Sample hyperparameters of Brownian Motion
        % tic;
        thetasamp = Logtheta(:,i);
        parfor j = 1:D % Use parfor to sample hyperparameters in parallel
			thetasamp(j) = SamplerTheta( UU(j,:), thetasamp(j), n, w, m, inv_BrownianMat ); % Sample hyperparameters by calling the slice sampler.
        end
        Logtheta(:,i+1) = thetasamp; % Save sampled hyperparameters.
        % toc 

       %% Sample copula model parameters using Spherical HMC
        % tic;
        [pmf1smdDATA,dpmfDATA] = pmfdpmfDATA( Y, exp(UU)./(1+exp(UU)) ); % Calculate the pmf and derivative of pmf
        [theta,acpyes] = SphHMC(theta, @(theta,der)Ut(theta,der,pmf1smdDATA,dpmfDATA), Stepsz, NLeap); % Sample theta by calling spherical HMC
        beta = abs(theta).*theta;
        beta = beta(1:K); % Map the sample on the sphere to the copula model parameters.
        Beta(:,i+1) = beta; % Save sampled copula model parameters
        accp = accp + acpyes; % Accumulate number of acceptance
        if i==BurnIn  
            disp('Burn in completed!');
        end                   	
        if i >= BurnIn % save sampled copula model parameters
        	Samp(:,i+1-BurnIn) = beta;
            acpt = acpt + acpyes; % Accumulate number of acceptance
            logwt(i+1-BurnIn) = K*log(2)+sum(log(abs(theta))); % log weights of samples
        end
        % toc
        
    end
    TIME = toc; % Time to end MCMC

end

        
        
