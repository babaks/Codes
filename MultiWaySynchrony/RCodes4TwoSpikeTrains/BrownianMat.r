#######################################################################################################################
## Function to calculate the covariance matrix of Brownian Motion
## Inputs: n is number of time bins
## Outputs: Covariance matrix of Brownian Motion
#######################################################################################################################
BrownianMat <- function(n){
	brownian.mat <- diag(1:n)
	for( i in 1:(n-1) ){
		brownian.mat[i,(i+1):n] <- rep( i, length((i+1):n) )
	}
	for( i in 1:(n-1) )
		for( j in (i+1):n )
			brownian.mat[j,i] <- brownian.mat[i,j]
	return( brownian.mat )	
}
#######################################################################################################################
## Function to calculate the Cholesky decomposition of covariance matrix of Brownian Motion
## Inputs: n is number of time bins
## Outputs: Cholesky decomposition of covariance matrix of Brownian Motion
#######################################################################################################################
BrownianChol <- function(n){
	brownian.chol <- diag(n)
	for( i in 1:(n-1) ){
		brownian.chol[i,(i+1):n] <- 1
	}
	return( brownian.chol )	
}