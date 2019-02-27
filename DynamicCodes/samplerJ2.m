%%
%%%%% Sampler to sample hyperparameters
%%

function newParam = samplerJ2( yy, logeta, logrho, logalpha, logJ2, n, w, m, DiffMat )
	% yy is the latent variable

	z = loglikeGpHyper( yy, logeta, logrho, logalpha, logJ2, n, DiffMat ) - exprnd(1);

	% Stepping out
	u = rand(1);
	L = logJ2 - w * u;
	R = L + w;
	v = rand( 1 );
	J = floor( m * v );
	K = ( m - 1 ) - J;

	while J>0 && z<loglikeGpHyper( yy, logeta, logrho, logalpha, L, n, DiffMat )
		L = L - w;
		J = J - 1;
    end
    
	while K>0 && z<loglikeGpHyper( yy, logeta, logrho, logalpha, R, n, DiffMat )
		R = R + w;
		K = K - 1;
    end
    
	% Shrinkage to obtain a new sample
	u = rand(1);
	newParam = L + u * ( R - L );

	while z > loglikeGpHyper( yy, logeta, logrho, logalpha, newParam, n, DiffMat )
        
		if newParam < logJ2
			L = newParam;
        else
			R = newParam;
        end
        
		u = rand(1);
		newParam = L + u*(R-L);
    end
end
