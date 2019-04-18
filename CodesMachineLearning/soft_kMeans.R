library(ggplot2)
library(mvtnorm)
library(bayesm)

# Generating some random data from 5 normal distributions 
n <- seq(60, 100, 10)
library(mvtnorm)
library(bayesm)
set.seed(5)
X <- rmvnorm(n[1], mean=rnorm(2, 0, 5), sigma=rwishart(3, diag(2))$IW)
for(i in 2:5){
  
  X <- rbind(X, rmvnorm(n[i], mean=rnorm(2, 0, 5), sigma=rwishart(5, diag(2))$IW))
  
}

N <- nrow(X)

data <- data.frame(x1 = X[, 1], x2 = X[, 2])
jpeg('../figures/softKmeansData.jpeg')
p <- ggplot(data, aes(x1, x2)) + xlim(c(-10, 15)) + ylim(-16, 8)
p + geom_point()
dev.off()


# Soft K-means clustering with K=5
K = 5

# Initializing the parameters 
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)
# I run the algorithm for 20 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:20){
# E-step   
for (k in 1:K){
  d[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])
}

total.d <- rowSums(d)  
for (k in 1:K){
  p[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])/total.d
}


# M-step
pi <- colSums(p)/N

sum.p <- colSums(p)  
for (k in 1:K){
  mu[k, ] <- colSums(p[, k]*X)/sum.p[k]
}

for (k in 1:K){
  temp <- 0
for(j in 1:N){
  temp <- temp + p[j, k]*(X[j, ] - mu[k, ])%*%t(X[j, ] - mu[k, ])
}
sigma[, , k] <- temp/sum.p[k]
}

  
  
}  


# Hard clustering: assigning each data point to the component with the highest probability
clus <- apply(p, 1, which.max)
data <- data.frame(x1 = X[, 1], x2 = X[, 2], clus = factor(clus))



# Simulating data from the fitted mixture model for comparison. 
Y <- NULL
for(k in 1:K){
  
  Y <- rbind(Y, rmvnorm(round(pi[k]*50000), mu[k, ], sigma[, , k]))
  
}

mix <- data.frame(y1 = Y[, 1], y2 = Y[, 2])


library(MASS)  
library(RColorBrewer)
k <- 10
my.cols <- rev(brewer.pal(k, "RdYlBu"))

z <- kde2d(Y[,1], Y[,2], n=100)


# # Plot using contour function
# pdf('../figures/softKmeansResults1.pdf')
# plot(X, xlab="X1", ylab="X2", pch=19, cex=.4)
# contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
# dev.off()

# Plot using ggplot2
jpeg('../figures/softKmeansResultsK5.jpeg')
pl <- ggplot(data, aes(x1, x2)) + xlim(c(-10, 15)) + ylim(-16, 8)
pl <- pl + geom_point(aes(color=clus, shape=clus, size=3)) 
pl <- pl + geom_density2d(aes(x=y1, y=y2), size=1, alpha=.3, data=mix ) + theme(legend.position="none")
pl
dev.off()






############



# Soft K-means clustering with K=3
K = 3

# Initializing the parameters 
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)
# I run the algorithm for 20 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:20){
# E-step  
for (k in 1:K){
  d[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])
}

total.d <- rowSums(d)  
for (k in 1:K){
  p[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])/total.d
}


# M-step
pi <- colSums(p)/N

sum.p <- colSums(p)  
for (k in 1:K){
  mu[k, ] <- colSums(p[, k]*X)/sum.p[k]
}

for (k in 1:K){
  temp <- 0
for(j in 1:N){
  temp <- temp + p[j, k]*(X[j, ] - mu[k, ])%*%t(X[j, ] - mu[k, ])
}
sigma[, , k] <- temp/sum.p[k]
}

  
  
}  


# Hard clustering: assigning each data point to the component with the highest probability
clus <- apply(p, 1, which.max)
data <- data.frame(x1 = X[, 1], x2 = X[, 2], clus = factor(clus))



# Simulating data from the fitted mixture model for comparison. 
Y <- NULL
for(k in 1:K){
  
  Y <- rbind(Y, rmvnorm(round(pi[k]*50000), mu[k, ], sigma[, , k]))
  
}

mix <- data.frame(y1 = Y[, 1], y2 = Y[, 2])


library(MASS)  

# # Plot using contour function

# library(RColorBrewer)
# k <- 10
# my.cols <- rev(brewer.pal(k, "RdYlBu"))
# 
# z <- kde2d(Y[,1], Y[,2], n=100)

# pdf('../figures/softKmeansResults1.pdf')
# plot(X, xlab="X1", ylab="X2", pch=19, cex=.4)
# contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
# dev.off()

# Plot using ggplot2
jpeg('../figures/softKmeansResultsK3.jpeg')
pl <- ggplot(data, aes(x1, x2)) + xlim(c(-10, 15)) + ylim(-16, 8)
pl <- pl + geom_point(aes(color=clus, shape=clus, size=3)) 
pl <- pl + geom_density2d(aes(x=y1, y=y2), size=1, alpha=.3, data=mix ) + theme(legend.position="none")
pl
dev.off()
