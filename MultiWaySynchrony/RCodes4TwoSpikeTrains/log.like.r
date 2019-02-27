#######################################################################################################################
## Log-Likelihood of P(t|y) (Binomial Distribution) where t is outcome and y is latent variable 
## This is log-likelihood for only one trial
## Inputs: t1 and t2 are the spike trains of the two neurons, p1 and p2 are the marginal firing probabilities of the
## two neurons, zeta is the exta term, lag is the lag between the two neurons, index is the index of trial, Bound is
## the Bound of lag.
## Outputs: Log-likelihood of one trial of spike trains.
#######################################################################################################################
loglike.single.trial <- function( t1, t2, p1, p2, zeta, lag, index, Bound ){

	#browser()
	nn <- length(p1)
	if( lag[index] == 0 ){ ## Two neurons are synchrony
		upper <- min( pmin( 1/p1, 1/p2 ) ) ## Upper bound of zeta
		lower <- max( pmax( (p1+p2-1)/(p1*p2), 0 ) ) ## Lower bound of zeta
		p11 <- p1 * p2 * ( (upper-lower)*zeta[0+Bound+1] + lower ) # Joint probability of cofiring
		p10 <- p1 - p11
		p01 <- p2 - p11
		p00 <- 1- p11 - p10 - p01
		if( sum(as.numeric(p00<=0))>0 )
			browser()

		log.like <- 0
		log.like <- log.like + sum( as.numeric( t1[,index]==1 & t2[,index]==1) * log(p11) )
		log.like <- log.like + sum( as.numeric( t1[,index]==1 & t2[,index]==0) * log(p10) )
		log.like <- log.like + sum( as.numeric( t1[,index]==0 & t2[,index]==1) * log(p01) )
		log.like <- log.like + sum( as.numeric( t1[,index]==0 & t2[,index]==0) * log(p00) )

	} else if( lag[index] > 0 ) ## Positive lags
	{
		p1.l <- p1[1:(nn-lag[index])]
		p2.l <- p2[(lag[index]+1):nn]
		upper <- min( pmin( 1/p1.l, 1/p2.l ) )
		lower <- max( pmax( (p1.l+p2.l-1)/(p1.l*p2.l), 0 ) )
		p11 <- p1.l * p2.l * ( (upper-lower)*zeta[lag[index]+Bound+1] + lower )
		p10 <- p1.l - p11
		p01 <- p2.l - p11
		p00 <- 1- p11 - p10 - p01
		if( sum(as.numeric(p00<=0))>0 )
			browser()

		log.like <- 0

		log.like <- log.like + sum( as.numeric( t1[1:(nn-lag[index]),index]==1 & t2[(1+lag[index]):nn,index]==1) * log(p11) )
		log.like <- log.like + sum( as.numeric( t1[1:(nn-lag[index]),index]==1 & t2[(1+lag[index]):nn,index]==0) * log(p10) )
		log.like <- log.like + sum( as.numeric( t1[1:(nn-lag[index]),index]==0 & t2[(1+lag[index]):nn,index]==1) * log(p01) )
		log.like <- log.like + sum( as.numeric( t1[1:(nn-lag[index]),index]==0 & t2[(1+lag[index]):nn,index]==0) * log(p00) )
		log.like <- log.like + sum( t1[(nn-lag[index]+1):nn,index]*log(p1[(nn-lag[index]+1):nn]) ) + sum( (1-t1[(nn-lag[index]+1):nn,index])*log(1-p1[(nn-lag[index]+1):nn]) )
		log.like <- log.like + sum( t2[1:lag[index],index]*log(p2[1:lag[index]]) ) + sum( (1-t2[1:lag[index],index])*log(1-p2[1:lag[index]]) )

	} else if( lag[index] < 0 ) ## Negative lags
	{
		p1.l <- p1[(1-lag[index]):nn]
		p2.l <- p2[1:(nn+lag[index])]
		upper <- min( pmin( 1/p1.l, 1/p2.l ) )
		lower <- max( pmax( (p1.l+p2.l-1)/(p1.l*p2.l), 0 ) )
		p11 <- p1.l * p2.l * ( (upper-lower)*zeta[lag[index]+Bound+1] + lower )
		p10 <- p1.l - p11
		p01 <- p2.l - p11
		p00 <- 1- p11 - p10 - p01
		if( sum(as.numeric(p00<=0))>0 )
			browser()

		log.like <- 0
		log.like <- log.like + sum( as.numeric( t1[(1-lag[index]):nn,index]==1 & t2[1:(nn+lag[index]),index]==1) * log(p11) )
		log.like <- log.like + sum( as.numeric( t1[(1-lag[index]):nn,index]==1 & t2[1:(nn+lag[index]),index]==0) * log(p10) )
		log.like <- log.like + sum( as.numeric( t1[(1-lag[index]):nn,index]==0 & t2[1:(nn+lag[index]),index]==1) * log(p01) )
		log.like <- log.like + sum( as.numeric( t1[(1-lag[index]):nn,index]==0 & t2[1:(nn+lag[index]),index]==0) * log(p00) )
		log.like <- log.like + sum( t1[1:(-lag[index]),index]*log(p1[1:(-lag[index])]) ) + sum( (1-t1[1:(-lag[index]),index])*log(1-p1[1:(-lag[index])]) )
		log.like <- log.like + sum( t2[(nn+lag[index]+1):nn,index]*log(p2[(nn+lag[index]+1):nn]) ) + sum( (1-t2[(nn+lag[index]+1):nn,index])*log(1-p2[(nn+lag[index]+1):nn]) )

	}

	return( log.like )

}

#######################################################################################################################
## Log-Likelihood of P(t|y) (Binomial Distribution) where t is outcome and y is latent variable 
## This is log-likelihood for two neurons (multiple trials)
## Inputs: t1 and t2 are the spike trains of the two neurons, p1 and p2 are the marginal firing probabilities of the
## two neurons, zeta is the exta term, lag is the lag between the two neurons, Bound is the Bound of lag.
## Outputs: Log-likelihood of two spike trains
#######################################################################################################################
loglike.data <- function( t1, t2, p1, p2, zeta, lag, Bound ){

	#browser()
	loglike <- 0
	for( g in 1:dim(t1)[2] )
		## Log-likelihood of one trial
		loglike <- loglike + loglike.single.trial( t1, t2, p1, p2, zeta, lag, index=g, Bound ) 
	return( loglike )

}

#######################################################################################################################
## Log-Likelihood w.r.t. latent variables (this is for one neuron)
## Inputs: y is the spike trains of one neuron, u is the latent variabls and invC is the inverse of covariance matrix
## of the Gaussian Process Prior.
## Outputs: Log-likelihood of data w.r.t. latent variables
#######################################################################################################################
loglike.latent <- function( u, y, invC ){ 

	#browser()
 	N <- dim(y)[2]
    	p = exp(u) / (1+exp(u)) ## Calculate firing probabilities
 	p.mat <- matrix(rep(p,N),ncol=N,byrow=FALSE)
   	loglik <- sum( y * log(p.mat) + (1-y) * log(1-p.mat) ) ## Calculate loglikelihood

	return( loglik )

}

#######################################################################################################################
## Log-likelihood of hyperparameters
## Inputs: y is the latent variables, log.theta is the log of the hyperparameter (theta), n is the number of time bins,
## inv.brownian.mat is the inverse of the covariance matrix of Brownian Motion.
## Outputs: Log-likelihood of hyperparameters
#######################################################################################################################
 
loglike.hyper <- function( y, log.theta, n, inv.brownian.mat ){

	
	log.prior <- -0.5*log.theta^2/3^2 ## Log-normal priors for hperparameters, log( theta^2 ) ~ N(0,3)
	 
	invC <- exp(log.theta) * inv.brownian.mat ## Inverse of Covariance Matrix
	logdetC <- -n*log.theta + log( det.brownian(1:n,n) ) ## Log of determinant of the inverse of covariance matrix
	log.like <- -0.5 * logdetC - 0.5 * t(y) %*% invC %*% y ## log-likelihood L=-n/2log(2pi)-1/2log(det(C))-1/2y^TC^(-1)y

	## log-posterior of hyperparameters
	log.post <- log.prior + log.like

	return( log.post )

}


