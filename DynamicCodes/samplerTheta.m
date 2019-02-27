%%
%%%%% Sampler to sample hyperparameter eta
%%

function newParam = samplerTheta( yy, logtheta, n, w, m, invBrownianMat )
	% yy is the latent variable

	z = loglikeBmHyper( yy, logtheta, n, invBrownianMat ) - exprnd(1);

	% Stepping out
	u = rand(1);
	L = logtheta - w * u;
	R = L + w;
	v = rand( 1 );
	J = floor( m * v );
	K = ( m - 1 ) - J;
    
    while J>0 && z<loglikeBmHyper( yy, L, n, invBrownianMat )
		L = L - w;
		J = J - 1;
    end
    while K>0 && z<loglikeBmHyper( yy, R, n, invBrownianMat )
		R = R + w;
		K = K - 1;
    end
    
	% Shrinkage to obtain a new sample
	u = rand(1);
	newParam = L + u * ( R - L );

    while z > loglikeBmHyper( yy, newParam, n, invBrownianMat )     
		
        if newParam < logtheta
			L = newParam;
        else
			R = newParam;
        end
        
		u = rand(1);
		newParam = L + u*(R-L);
    end
end

        
        
