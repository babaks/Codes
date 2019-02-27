%%
%%%%% Sampler to sample hyperparameter eta
%%

function newParam = samplerEta( yy, logeta, logrho, logalpha, logJ2, n, w, m, DiffMat )
	% yy is the latent variable

	z = loglikeGpHyper( yy, logeta, logrho, logalpha, logJ2, n, DiffMat ) - exprnd(1);

	% Stepping out
	u = rand(1);
	L = logeta - w * u;
	R = L + w;
	v = rand( 1 );
	J = floor( m * v );
	K = ( m - 1 ) - J;

	while J>0 && z<loglikeGpHyper( yy, L, logrho, logalpha, logJ2, n, DiffMat )
		L = L - w;
		J = J - 1;
    end
    
	while K>0 && z<loglikeGpHyper( yy, R, logrho, logalpha, logJ2, n, DiffMat )
		R = R + w;
		K = K - 1;
    end
    
	% Shrinkage to obtain a new sample
	u = rand(1);
	newParam = L + u * ( R - L );

	while z > loglikeGpHyper( yy, newParam, logrho, logalpha, logJ2, n, DiffMat )
        
		if newParam < logeta
			L = newParam;
        else
			R = newParam;
        end
        
		u = rand(1);
		newParam = L + u*(R-L);
    end
end

        
        
