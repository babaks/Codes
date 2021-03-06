%%
%%%%% Sampler to sample hyperparameters rho
%%

function newParam = samplerRho( yy, logeta, logrho, logalpha, logJ2, n, w, m, DiffMat )
	% yy is the latent variable

	z = loglikeGpHyper( yy, logeta, logrho, logalpha, logJ2, n, DiffMat ) - exprnd(1);

	% Stepping out
	u = rand(1);
	L = logrho - w * u;
	R = L + w;
	v = rand( 1 );
	J = floor( m * v );
	K = ( m - 1 ) - J;

	while J>0 && z<loglikeGpHyper( yy, logeta, L, logalpha, logJ2, n, DiffMat )
		L = L - w;
		J = J - 1;
    end
    
	while K>0 && z<loglikeGpHyper( yy, logeta, R, logalpha, logJ2, n, DiffMat )
		R = R + w;
		K = K - 1;
    end
    
	% Shrinkage to obtain a new sample
	u = rand(1);
	newParam = L + u * ( R - L );

	while z > loglikeGpHyper( yy, logeta, newParam, logalpha, logJ2, n, DiffMat )
        
		if newParam < logrho
			L = newParam;
        else
			R = newParam;
        end
        
		u = rand(1);
		newParam = L + u*(R-L);
    end
end
