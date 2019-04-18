library(ggplot2)

# We generate data from two normal distributions: 20 samples with means 0 and 10 samples with 5; we fist the variance at 1
set.seed(1)
x <- c(rnorm(20, 0, 1), rnorm(10, 5, 1))
n = length(x)

# We generate 10000 MCMC samples
nIter = 10000

# We set gamma to 3. In genral, however, it's better to put a prior on gamma
gamma <- 3

# Baeline distribution G0 is assumed to be N(0, 10)
m = 0
tau2 <- 10

# Here, we fist the variance at 1 so we are not sampling it. 
sigma2 <- 1

# Theta contains the means for each observation; note that these means are not unique. 
theta <- matrix(0, n, nIter)


for(j in 1:nIter){
  for(i in 1:n){
    # Probability (up to some constant) that the mean comes from G0
    q0 <- gamma*dnorm(x[i], m, sigma2+tau2) 
    
    # Probability that the mean for the ith observation is the same as the mean of the jth observation
    qj <- dnorm(x[i], theta[-i, j], 1) 
    
    # Normalized probabilities
    p = c(q0, qj)/sum(c(q0, qj)) 
    
    # We sample a new theta from the posterior distribution given G0 and create a vector of possible means to choose from; the elements of this vector correspond to p
    Thetas <- c( rnorm(1, (m/tau2+x[i]/sigma2)/(1/tau2+1/sigma2), 1/(1/tau2+1/sigma2) ), theta[-i, j]) 
    
    # we then sample a theta with probabilty p
    theta[i, j] <- sample(Thetas, 1, prob=p) 
  }
}



# We throguh away the first 3000 samples and take average using the remaining samples
plot(rowMeans(theta[, 3001:nIter]))


# Simulating data using the MCMC samples
Y = NULL
count = 0
for(i in seq(5001, 10000, 500)){
  count = count + 1
  Y = rbind(Y, cbind(rnorm(10*n, theta[, i], 1), count))
  
}


# Ploting 10 posterior samples
jpeg('../figures/dpUniDensity.jpeg')
simData.df <- data.frame(x=Y[, 1], y = Y[, 2])
origData.df <- data.frame(x = x, y=rep(0, n))
p <- ggplot(simData.df, aes(x=x, group=y)) + geom_density() # +theme(legend.position="none")
p <- p + geom_point(aes(x=x, y=y), size=4, cex=16, data=df)
p
dev.off()



# This is a highly vectorized version of the algorithm
# for(j in 1:nIter){
#   for(i in 1:n){   
#     p = c(q0 <- gamma*dnorm(y[i], m, sigma2+tau2), qj <- dnorm(y[i], theta[-i, j], 1))/sum(c(q0, qj))
#     Thetas <- c( rnorm(1, (m/tau2+y[i]/sigma2)/(1/tau2+1/sigma2), 1/(1/tau2+1/sigma2) ), theta[-i, j])
#     theta[i, j] <- sample(Thetas, 1, prob=p)
#   }
# }
