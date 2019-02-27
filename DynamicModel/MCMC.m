%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to do mcmc sampling for latent variables, hyperparameters of Brownian Motion, and parameters of the copula model using elliptical slice 
% sampler, slice sampler, and so forth.
% Input: Y is the data (Spike Trains), Niteration is the number of iterations of MCMC, w and m are step size and step number of slice sampler.
% Output: Logtheta, Logrho, Logalpha and LogJ2 are the log of the hyperparameter of Gaussian Process priors.
% TR  is the sampled thresholds and II is the sampled indicators.  
% logwt is the log weights of the spherical HMC samples, acpt is the acceptance rate of spherical HMC, and TIME is the time to do MCMC sampling.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Logeta, Logrho, Logalpha, LogJ2, TR, UU_samp, II, SIGMA, TIME] = MCMC( time, Y, Niteration, w, m  )	
    [D,n,N] = size(Y); % dimension of data; N: number of trials; n: number of time bins; D: number of neurons;
    
    %% Initialization for hyperparameters and latent variables
	Logeta = zeros( D, Niteration );
	Logrho = zeros( D, Niteration );
	Logalpha = zeros( D, Niteration );
	LogJ2 = zeros( D, Niteration );
    TR = zeros( n, D, Niteration ); 
    Tr_current = zeros(n,D); 
    UU_samp = zeros( D, n, N, Niteration );
    SIGMA = zeros( D, D, Niteration );
    for i = 1:D 
        SIGMA(i,i,1) = 1;
    end    
    SIGMA_current = zeros(D,D);
    SIGMA_current(:,:) = SIGMA(:,:,1);
    Chol_Sigma = chol(SIGMA_current);
    shrinkage = -1;
    pt = zeros(D,n,Niteration)+2*normcdf(exp(shrinkage))-1;
	pt_current = zeros(D,n);
    qt_current = zeros(D,n) + shrinkage; 
    II = zeros( D, n, Niteration ); 
    II_current = zeros(D,n);
    Logtheta = zeros(D,Niteration) + 3;
 
    DiffMat = zeros(n,n); % used to calculate the covariance matrix for Gaussian Process
    for i = 1:n
        for j = 1:n
            DiffMat(i,j)=(time(i)-time(j))^2;
        end
    end
	inv_BrownianMat = SolveBrownian(1:n,n); %% Inverse of covariance matrix of Brownian Mothion.
    brownian_Chol = BrownianChol(n);   
    
    for i = 1:(Niteration-1)
        
        disp( i );     %% Index of iteration

        %% Sample latent variables
        tic;
        II_current(:,:) = II(:,:,i);
        Tr_current(:,:) = TR(:,:,i);
        SIGMA_current(:,:) = SIGMA(:,:,i);
        UU = samplerLatent( Y, SIGMA_current, Tr_current, II_current );
        UU_samp(:,:,:,i) = UU;
        toc
        
        %% Sample Thresholds
        tic;
        parfor j=1:D % Use parfor to sample latent variables in parallel
            minU = zeros(n,1);
            maxU = zeros(n,1);
            utemp = zeros(n,N);
            utemp(:,:) = UU(j,:,:); % data of the j^th neuron
            ytemp = zeros(n,N);
            ytemp(:,:) = Y(j,:,:);
            for k = 1:n
                %index1 = find( ytemp(k,:)==1 );
                %index0 = find( ytemp(k,:)==0 );
                index1 = ytemp(k,:)==1;
                index0 = ytemp(k,:)==0;
                if sum(index1) > 0
                    minU(k) = max(utemp(k,index1));
                else
                    minU(k) = -Inf;
                end
                if sum(index0) > 0
                    maxU(k) = min(utemp(k,index0));
                else
                    maxU(k) = Inf;
                end                
            end
			C = exp(LogJ2(j,i)) * diag(ones(1,n)) + exp(Logeta(j,i)) * exp( -DiffMat * exp(Logrho(j,i)) ) + exp(Logalpha(j,i))*ones(n);            
			cholC = chol(C); %% Cholesky decomposition
            Tr_current(:,j) = samplerTr( Tr_current(:,j), minU, maxU, n, cholC ); % Sample latent variables by calling the elliptical slice sampler.
        end
        TR(:,:,i+1) = Tr_current;        
        toc
       
        %% Sampler GP Hyperparameters       
        tic;
        etasamp = Logeta(:,i);
        rhosamp = Logrho(:,i);
        alphasamp = Logalpha(:,i);
        J2samp = LogJ2(:,i);
        parfor j = 1:D
           %% sample hyper-parameters
			etasamp(j) = samplerEta( Tr_current(:,j), etasamp(j), rhosamp(j), alphasamp(j), J2samp(j), n, w, m, DiffMat );
			rhosamp(j) = samplerRho( Tr_current(:,j), etasamp(j), rhosamp(j), alphasamp(j), J2samp(j), n, w, m, DiffMat );
			alphasamp(j) = samplerAlpha( Tr_current(:,j), etasamp(j), rhosamp(j), alphasamp(j), J2samp(j), n, w, m, DiffMat );
			J2samp(j) = samplerJ2( Tr_current(:,j), etasamp(j), rhosamp(j), alphasamp(j), J2samp(j), n, w, m, DiffMat );
        end
        Logeta(:,i+1) = etasamp;
        Logrho(:,i+1) = rhosamp;
        Logalpha(:,i+1) = alphasamp;
        LogJ2(:,i+1) = J2samp;
        toc
        
        %% Sample Correlation Matrix
        tic;
        Chol_Sigma = samplerSigma( UU, Chol_Sigma, II_current, D ); 
        toc
        SIGMA(:,:,i+1) = Chol_Sigma' * Chol_Sigma;
        
        %% Sample Indicators
        tic;
    	SIGMA_current(:,:) = SIGMA(:,:,i+1);
        pt_current(:,:) = pt(:,:,i);
        index_pair = 1:D;
        for j = 1:(D/2) 
			[output1,output2] = samplerInd ( UU, SIGMA_current, index_pair(2*j-1), index_pair(2*j), II_current, pt_current ); 
            II_current(index_pair(2*j-1),:) = output1;
            II_current(index_pair(2*j),:) = output2;
        end
        II(:,:,i+1) = II_current; % Save sampled indicators.
        toc
        
        %% Sample pt
        %% Sample probabilities of being involved in the network and BM
        %% hyperparameters
        tic;
        thetasamp = Logtheta(:,i);
        parfor j =2:D
            cholC = sqrt(1/exp(thetasamp(j))) * brownian_Chol;
            qt_current(j,:) = samplerQt( qt_current(j,:), II_current(j,:), n, cholC, shrinkage );
            thetasamp(j) = samplerTheta( qt_current(j,:)-shrinkage, thetasamp(j), n, w, m, inv_BrownianMat );
        end
        pt(:,:,i+1) = 2*normcdf(exp(qt_current))-1; 
        Logtheta(:,i+1) = thetasamp;  
        toc
        
    end
    
    TIME = toc; % Time to end MCMC

end

        
        
