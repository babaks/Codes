#######################################################################################################################
## Function to sample latent variables using elliptical slice sampler
## Inputs: u.current is the current state of latent variables, y is the spike train data, n is number of time bins, 
## chol.C is the Cholesky decompostion of the covariance matrix of Gaussian Process prior and invC is the inverse 
## of the covariance matrix of Gaussian Process prior.
## Outputs: New state of hyperparameters
#######################################################################################################################
sampler.latent <- function( u.current, y, n, chol.C, invC ){
	
	# browser()
	v <- t(chol.C)  %*% rnorm(n,0,1) ## Propose a new state of latent variables
	u <- runif( 1, 0, 1 )
	logy <- loglike.latent( u.current, y, invC ) + log(u) ## Calculate the threshold for elliptical slice sampler
	
	theta <- runif( 1, 0, 2*pi ) ## Randomly generate an angle
	theta.min <- theta - 2*pi ## Define a bracket
	theta.max <- theta
	
	u.prop <- u.current * cos(theta) + v * sin(theta) ## Proposed state of latent variables
	while ( logy > loglike.latent( u.prop, y, invC ) ){
		if ( theta < 0 ) { ## Shrinkage
			theta.min <- theta
		}else{
			theta.max <- theta
		}
		theta <- runif( 1, theta.min, theta.max ) ## Generate a new angle
		u.prop <- u.current * cos(theta) + v * sin(theta) ## Propose a new state for latent variables
	}
	return( u.prop )
}

#######################################################################################################################
## Function to sample the extra term zeta using Metropolis Hastings sampling.
## Inputs: U and V are the latent variables for two neurons, W.before is the current state, lag is the lag between 
## the two neurons, t1 and t2 are the spike trains of the two neurons, epsilon is the step size and Bound is the 
## bound of lags. 
## Outputs: New state of zeta
#######################################################################################################################
sampler.w <- function( U, V, W.before, lag, t1, t2, epsilon, Bound ){
	
	p1 <- exp( U) / ( 1+exp(U) ) ## Calculate marginal firing probabilities
	p2 <- exp(V) / ( 1+exp(V) )
	zeta <- exp(W.before) / ( 1+exp(W.before) ) 
	## Calculate log-likelihood before sampling
	loglike.before.1 <- loglike.data( t1, t2, p1, p2, zeta, lag ) ## Loglikelihood of the current state
	## Sampler using M-H 	
	W.prop <- W.before + epsilon * rnorm( 1 )
	zeta.prop <- exp( W.prop) / ( 1+exp(W.prop) )
	loglike.prop.1 <- loglike.data( t1, t2, p1, p2, zeta.prop, lag ) ## Loglikelihood of the propose state

	## Calculate acceptance probability
	a <- min( 1, exp( loglike.prop.1 - loglike.before.1 ) ) ## The proposal is symmetric
		
	## Decide whether the proposal is accepted or not 
	u <- runif(1,0,1)
	W.sample <- W.before
	if( u < a ){
		 W.sample <- W.prop
	}
	return( W.sample )

}


#######################################################################################################################
## Function to sample the lag using Metropolis Hastings sampling.
## Inputs: U and V are the latent variables for two neurons, W is the extra term zeta, lag is the current state of 
## lag between the two neurons, t1 and t2 are the spike trains of the two neurons, lag.hyper is the prior probabilities
## of lages, and Bound is the bound of lags. 
## Outputs: New state of lag
#######################################################################################################################

sampler.lag <- function( U, V, W, t1, t2, lag, lag.hyper, Bound ){

	p1 <- exp( U) / ( 1+exp(U) ) ## Marginal firing probabilties
	p2 <- exp(V) / ( 1+exp(V) )
	zeta <- exp(W) / ( 1+exp(W) ) ## Extra term zeta

	lag.sample <- rep( 0, length(lag) )
	for( s in 1:dim(t1)[2] ){
		p.lag <- numeric(0)		
		for( E in -Bound:Bound){

			lag[s] <- E
			p.lag <- c( p.lag, loglike.single.trial( t1, t2, p1, p2, zeta, lag, index=s, Bound ) )
		}
	
		p.lag <- exp(p.lag-max(p.lag))/sum(exp(p.lag-max(p.lag)))
		Prob <- p.lag*lag.hyper/sum(p.lag*lag.hyper)

		lag.samp <- rmultinom( 1, size=1, prob=Prob ) ## Sample lag from multinomial distribution
		lag[s] <- lag.sample[s] <- which(lag.samp==1) - (Bound+1)
	}
	
	return( lag.sample )
}


#######################################################################################################################
## Function to sample prior probabilties of lags.
## Inputs: lag is the current state of lag between the two neurons, and Bound is the bound of lags. 
## Outputs: New state of prior probabilties of lags
#######################################################################################################################
sampler.lag.hyper <- function( lag, Bound ){

	X <- rep( 0, 2*Bound+1 )
	for( s in -Bound:Bound )
		X[s+Bound+1] <- sum( as.numeric(lag==s) )

	samp.lag.hyper <- as.numeric( rDirichlet.acomp( 1, X+0.01 ) ) ## Sample prior probabilties from Dirichlet distribution
	return( samp.lag.hyper )
}


