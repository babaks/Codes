#######################################################################################################################
## Function to sample hyperparameters using slice sampler
## Inputs: y is the latent variables, log.theta is the current state of hyperparameters, n is number of time bins, 
## inv.brownian.mat is the inverse of the covariance matrix of Brownian Motion, w and m are the step size and step 
## number of slikce sampler
## Outputs: New state of hyperparameters
#######################################################################################################################
sampler.theta <- function( y, log.theta, n, inv.brownian.mat, w, m ){

	z <- loglike.hyper( y, log.theta, n, inv.brownian.mat ) - rexp( 1 ) ## Threshold for slice sampler

	## Stepping out
	u <- runif(1)
	L <- log.theta - w * u
	R <- L + w
	v <- runif( 1 )
	J <- floor( m * v )
	K <- ( m - 1 ) - J

	while ( J>0 && z < loglike.hyper( y, L, n, inv.brownian.mat ) ){
		L <- L - w	
		J <- J - 1
	}
	while ( K>0 && z < loglike.hyper( y, R, n, inv.brownian.mat ) ){
		R <- R + w
		K <- K - 1
	}

	## Shrinkage to obtain a new sample
	u <- runif(1)
	newParam <- L + u * ( R - L )

	while ( z > loglike.hyper( y, newParam, n, inv.brownian.mat ) ){
		if ( newParam < log.theta ) {
			L <- newParam
		}else{
			R <- newParam
		}
		u <- runif(1)
		newParam <- L + u*(R-L)
	}
	
	return( newParam )
}

