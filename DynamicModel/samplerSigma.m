function output = samplerSigma( uu, Chol_bef, Ind, D )
%%
%% Function to draw samples for the covariance matrix using elliptical
%% slice sampler.
%%

v = zeros(D,D);
for ii = 1:D 
    v(ii,ii) = normrnd(0,3,1); % Propose a state from the exp prior
end
u = rand(1);
Sigma_bef = Chol_bef' * Chol_bef;
logy = loglikeSigma( uu, Sigma_bef, Ind ) + log(u); % Calculate the threshold for slice sampler

theta = unifrnd( 0, 2*pi, 1 ); % Draw a proposal
theta_min = theta - 2*pi; % Define a bracket
theta_max = theta; % Define a bracket

Chol_prop = Chol_bef;
for ii =1:D
    Chol_prop(ii,ii) = exp( log(Chol_bef(ii,ii)) * cos(theta) + v(ii,ii) * sin(theta) );  % Proposed a new state
end 
Sigma_prop = Chol_prop' * Chol_prop;
while logy > loglikeSigma( uu, Sigma_prop, Ind )% Compare the log-likelihood of the proposed state to the threshold
	
    if theta < 0 % Shrinkage
        theta_min = theta; 
    else
        theta_max = theta;
    end
    theta = unifrnd( theta_min, theta_max, 1); % Propose a new angle
    for ii =1:D
        Chol_prop(ii,ii) = exp( log(Chol_bef(ii,ii)) * cos(theta) + v(ii,ii) * sin(theta) );  % Proposed a new state
    end
    Sigma_prop = Chol_prop' * Chol_prop;
end

v = zeros(D,D);
for ii = 1:(D-1)
    for jj = (ii+1):D
        v(ii,jj) = normrnd(0,3,1); % Propose a state from the Laplace prior
    end
end
u = rand(1);
Chol_bef = Chol_prop;
Sigma_bef = Chol_bef' * Chol_bef;
logy = loglikeSigma( uu, Sigma_bef, Ind ) + log(u); % Calculate the threshold for slice sampler

theta = unifrnd( 0, 2*pi, 1 ); % Draw a proposal
theta_min = theta - 2*pi; % Define a bracket
theta_max = theta; % Define a bracket

Chol_prop = Chol_bef;
for ii = 1:(D-1)
    for jj = (ii+1):D
        Chol_prop(ii,jj) = Chol_bef(ii,jj) * cos(theta) + v(ii,jj) * sin(theta); % Propose a state
    end
end
Sigma_prop = Chol_prop' * Chol_prop;
while logy > loglikeSigma( uu, Sigma_prop, Ind )% Compare the log-likelihood of the proposed state to the threshold
	
    if theta < 0 % Shrinkage
        theta_min = theta; 
    else
        theta_max = theta;
    end
    theta = unifrnd( theta_min, theta_max, 1); % Propose a new angle
    for ii = 1:(D-1)
        for jj = (ii+1):D
            Chol_prop(ii,jj) = Chol_bef(ii,jj) * cos(theta) + v(ii,jj) * sin(theta); % Propose a state
        end
    end
    Sigma_prop = Chol_prop' * Chol_prop;
end

output = Chol_prop;   

end