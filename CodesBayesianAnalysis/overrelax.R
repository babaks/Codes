overrelaxGibbs <- function(rho= 0.99, alpha = -0.9, nIter=1000){
	
	x1 = rep(NA, nIter)
	x2 = rep(NA, nIter)
	
	x1[1] = 1
	x2[1] = 1


	x1.or = rep(NA, nIter)
	x2.or = rep(NA, nIter)
	
	x1.or[1] = 1
	x2.or[1] = 1
	
	for(i in 1:nIter){
		x1[i+1] = rnorm(1, rho*x2[i], sqrt(1-rho^2))
		x2[i+1] = rnorm(1, rho*x1[i+1], sqrt(1-rho^2))
		x1.or[i+1] = rho*x2.or[i] + alpha*(x1.or[i] - rho*x2.or[i])+sqrt(1-rho^2)*sqrt(1-alpha^2)*rnorm(1)
		x2.or[i+1] = rho*x1.or[i+1] + alpha*(x2.or[i] - rho*x1.or[i+1])+sqrt(1-rho^2)*sqrt(1-alpha^2)*rnorm(1)
		}
	
	return(list(x1=x1, x2=x2, x1.or = x1.or, x2.or=x2.or))
	
	}
	
	
res = overrelaxGibbs(nIter=5000)	


par(mfrow=c(2, 1))

plot(res$x1, type='l', ylab=expression(x[1]), xlab='Iteration', main='Without overrelaxation')

plot(res$x1.or, type='l', ylab=expression(x[1]), xlab='Iteration', main='With overrelaxation')



library(mvtnorm)

mu = c(0, 0)
sig = diag(2)
sig[1, 2] = 0.99
sig[2, 1] = 0.99
x=seq(-3, 3, .2)
y=seq(-3, 3, .2)

model <- function(a, b){dmvnorm(cbind(a, b), mu, sig)}
z = outer(x, y, model)

par(mfrow=c(2, 1))
contour(x, y, z, levels=0.01, main='Without overrelaxation', lty=2, col='darkgray', xlab=expression(x[1]), ylab=expression(x[2]))
lines(res$x1[1:100], res$x2[1:100])

contour(x, y, z, levels=0.01, main='With overrelaxation', lty=2, col='darkgray', xlab=expression(x[1]), ylab=expression(x[2]))
lines(res$x1.or[1:100], res$x2.or[1:100])
