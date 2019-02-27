%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to sample latent variables using elliptical slice sampler
% Inputs: u_current is the current state of latent variables, y is the data, n is number of time bins and cholC is the cholesky decomposition of covariance matrix
% Output: u_prop is the proposed state of latent variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_prop = SamplerLatent( u_current, y, n, cholC )
	
    v = normrnd(0,1,[1,n]) * cholC; % Propose a state from the GP prior
	u = rand(1);
	logy = LoglikeLatent( u_current, y ) + log(u); % Calculate the threshold for slice sampler
    
    theta = unifrnd( 0, 2*pi, 1 ); % Draw a proposal
	theta_min = theta - 2*pi; % Define a bracket
	theta_max = theta; % Define a bracket

	u_prop = u_current * cos(theta) + v * sin(theta); % Proposed state of latent variables
    while logy > LoglikeLatent( u_prop, y ) % Compare the log-likelihood of the proposed state to the threshold
	
        if theta < 0 % Shrinkage
			theta_min = theta; 
        else
			theta_max = theta;
        end
        
		theta = unifrnd( theta_min, theta_max, 1); % Propose a new angle
        u_prop = u_current * cos(theta) + v * sin(theta); % Propose a new state of latent variables
    end
end