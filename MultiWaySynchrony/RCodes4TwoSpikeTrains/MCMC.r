
#################################################################################################################################
## Function to do MCMC sampling by calling sub-functions
## Inputs: x is the time, t1 and t2 are the spike trains of the two neurons, Niteration is number of iterations, w and
## m is the step size and step length of slice sampler, epsilon is the step size of Metropolis Hastings, and Bound is 
## the bound of lags. 
## Outputs: samples of hyperparameters, samples of marginal firing probabilities (derived from latent variables), 
## samples of the extra term zeta, samples of lags, samples of prior probabilities of lags, and time to do MCMC sampling. 
#################################################################################################################################
BM_LAG_TTRE <- function( x, t1, t2, Niteration, w, m, epsilon, Bound=20 ){
	
	## Initialization for hyperparameters, latent variables, extra term and lags.
	n <- length( x )
	logtheta.u <- rep( 5, Niteration )
	logtheta.v <- rep( 5, Niteration )
	UU <- rnorm( n, -4, 0.1 )
	VV <- rnorm( n, -4, 0.1 )
	WW <- rnorm( 2*Bound+1, 0, 0.1 )
	P.u <- numeric(0) ## to store the predictive probability ( exp(y)/1+exp(y) )
	P.v <- numeric(0)
	Lag <- matrix( 0, nrow=dim(t1)[2], ncol=Niteration )
	Lag.hyper <- matrix( 1/(2*Bound+1), nrow=2*Bound+1, ncol=Niteration )
	ZETA <- numeric(0)

	brownian.chol <- BrownianChol(n) ## Cholesky decomposition of the covariance matrix of Brownian Motion		
	inv.brownian.mat <- solve.brownian( 1:n, n ) ## Inverse of the covariance matrix of Brownian Motion

	start <- proc.time()
	for( i in 1:(Niteration-1) ){

		#browser() 
		## Display index of iteration, and samples of hyperparameters 
		print( i )
		print( c(logtheta.u[i], UU[1:5]) )
		print( c(logtheta.v[i], VV[1:5]) )
		
		## sample latent variables 
		chol.C.u <- sqrt(1/exp(logtheta.u[i])) * brownian.chol ## Cholesky decomposition of the covariance matrix of Brownian Motion
		invC.u <- exp(logtheta.u[i]) * inv.brownian.mat ## Inverse of Covariance Matrix
		UU <- sampler.latent( UU, t1, n, chol.C.u, invC.u ) ## Sample latent variables by calling sub-functions
		p.u <- exp(UU)/(1+exp(UU)) ## Store the firing probabilities
		P.u <- cbind( P.u, p.u )
		
		chol.C.v <- sqrt(1/exp(logtheta.v[i])) * brownian.chol ## Cholesky decomposition of the covariance matrix of Brownian Motion
		invC.v <- exp(logtheta.v[i]) * inv.brownian.mat ## Inverse of Covariance Matrix
		VV <- sampler.latent( VV, t2, n, chol.C.v, invC.v ) ## Sample latent variables by calling sub-function
		p.v <- exp(VV)/(1+exp(VV)) ## Store the firing probabilities
		P.v <- cbind( P.v, p.v )

		## sample the extra term zeta
		WW.samp <- sampler.w( UU, VV, WW, Lag[,i], t1, t2, epsilon, Bound ) ## Sample extra term zeta by calling sub-function
		Upper <- min( pmin( 1/p.u, 1/p.v ) )
		Lower <- max( pmax( (p.u+p.v-1)/(p.u*p.v), 0 ) )
		ZETA <- c( ZETA, exp(WW)/(1+exp(WW))*(Upper-Lower) + Lower ) ## Store extra term zeta

		## sample Lag
		Lag[,i+1] <- sampler.lag( UU, VV, WW, t1, t2, Lag[,i], Lag.hyper[,i], Bound ) ## Sample lag by calling sub-function
		Lag.hyper[,i+1] <- sampler.lag.hyper( Lag[,i+1], Bound ) ## Store samples of lags

		## sample hyper-parameters
		logtheta.u[i+1] <- sampler.theta( UU, log.theta=logtheta.u[i], n, inv.brownian.mat, w, m ) ## Sample hyperparameters by calling sub-functions
		logtheta.v[i+1] <- sampler.theta( VV, log.theta=logtheta.v[i], n, inv.brownian.mat, w, m ) ## Sample hyperparameters by calling sub-functions	
	}

	Time = proc.time()-start ## Time to do MCMC sampling
	print( Time[1] )

	out <- list( logtheta.u, logtheta.v, P.u, P.v, ZETA, Lag, Lag.hyper, Time  )
	names( out ) <- c( "logtheta.u", "logtheta.v", "P1", "P2", "ZETA", "Lag", "Lag.hyper", "Time" )
	return( out )
}