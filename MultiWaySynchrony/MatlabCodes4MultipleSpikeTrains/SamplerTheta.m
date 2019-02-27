%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to sample hyperparameter using slice sampler
% Inputs: yy is the latent variables, logtheta is the current state of log of hyperparameter, n is number of time bins, w and m is the step size and step number, inv_BrownianMat is the inverse of covariance matrix
% Output: Proposed state of hyperparameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function newParam = SamplerTheta( yy, logtheta, n, w, m, inv_BrownianMat )

	z = LoglikeHyper( yy, logtheta, n, inv_BrownianMat ) - exprnd(1); % Thereshold for slice sampler

	% Stepping out
	u = rand(1);
	L = logtheta - w * u;
	R = L + w;
	v = rand( 1 );
	J = floor( m * v );
	K = ( m - 1 ) - J;    
    while J>0 && z<LoglikeHyper( yy, L, n, inv_BrownianMat )
		L = L - w;
		J = J - 1;
    end
    while K>0 && z<LoglikeHyper( yy, R, n, inv_BrownianMat )
		R = R + w;
		K = K - 1;
    end
    
	% Shrinkage to obtain a new sample
	u = rand(1);
	newParam = L + u * ( R - L );

    while z > LoglikeHyper( yy, newParam, n, inv_BrownianMat )     
		
        if newParam < logtheta
			L = newParam;
        else
			R = newParam;
        end
        
		u = rand(1);
		newParam = L + u*(R-L);
    end
end

        
        
