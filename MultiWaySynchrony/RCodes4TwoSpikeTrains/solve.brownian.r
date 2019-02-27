#######################################################################################################################
## Function to calculate the inverse of the covariance matrix of Brownian Motion
## Inputs: t is the time and n is number of time bins
## Outputs: inverse of the covariance matrix of Brownian Motion
#######################################################################################################################
solve.brownian <- function( t, n ){

	inv.Sigma <- matrix( 0, nrow=n, ncol=n )
	diag(inv.Sigma) <- c( 1/(t[2:n]-t[1:(n-1)])+1/(t[1:(n-1)]-c(0,t[1:(n-2)])), 1/(t[n]-t[n-1]) )
	inv.Sigma[1:(n-1),2:n] <- inv.Sigma[1:(n-1),2:n] - diag(1/(t[2:n]-t[1:(n-1)]))
	inv.Sigma[2:n,1:(n-1)] <- inv.Sigma[2:n,1:(n-1)] - diag(1/(t[2:n]-t[1:(n-1)]))
	
	return( inv.Sigma )
}

#######################################################################################################################
## Function to calculate the determinant of the covariance matrix of Brownian Motion
## Inputs: t is time and n is number of time bins
#######################################################################################################################
det.brownian <- function( t, n ){
	
	det.Sigma <- exp( sum(log(t[1:n]-c(0,t[1:(n-1)]))) )
	return( det.Sigma )
}


