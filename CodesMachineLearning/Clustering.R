
# Generating some random data from 5 normal distributions 
n <- seq(60, 100, 10)
library(mvtnorm)
library(bayesm)
set.seed(111)
X <- rmvnorm(n[1], mean=rnorm(2, 0, 5), sigma=rwishart(3, diag(2))$IW)
for(i in 2:5){
  
  X <- rbind(X, rmvnorm(n[i], mean=rnorm(2, 0, 5), sigma=rwishart(5, diag(2))$IW))
  
}

N <- nrow(X)

data <- data.frame(x1 = X[, 1], x2 = X[, 2])
#jpeg('../figures/kmeansData.jpeg')
p <- ggplot(data, aes(x1, x2))
p + geom_point()
#dev.off()


# K-means clustering with K=5
K = 5

# Initializing the centroids 
Cent <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
D <- matrix(NA, N, K)

# I run the algorithm for 100 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:100){
# Calculating the distance of each data point from the K centroids.   
for (k in 1:K){
  D[, k] <- rowSums( (sweep(X, 2, Cent[k, ]))^2)
}

# Assigning the data points to the closest centroid  
clus <- apply(D, 1, which.min)

# Finding the new centroids   
Cent <- by(X, INDICES=clus, FUN=colMeans)  
Cent <- do.call(rbind, Cent)  

}  


data <- data.frame(x1 = X[, 1], x2 = X[, 2], clus = factor(clus))


#jpeg('../figures/kmeansResults.jpeg')
p <- ggplot(data, aes(x1, x2))
p + geom_point(aes(color=clus, shape=clus)) 
#dev.off()




######## K-means clustering with K = 3

K = 3
Cent <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
D <- matrix(NA, N, K)

for(i in 1:10){
for (k in 1:K){
  D[, k] <- rowSums( (sweep(X, 2, Cent[k, ]))^2)
}

clus <- apply(D, 1, which.min)

Cent <- by(X, INDICES=clus, FUN=colMeans)  
Cent <- do.call(rbind, Cent)  

}  


data <- data.frame(x1 = X[, 1], x2 = X[, 2], clus = factor(clus))


#jpeg('../figures/kmeansResults2.jpeg')
p <- ggplot(data, aes(x1, x2))
p + geom_point(aes(color=clus, shape=clus)) 
#dev.off()





##### Hierarchical clustering

hc <- hclust(dist(X), "single")
dend <- as.dendrogram(hc)

plot(dend, leaflab = "none")
abline(h=4, lty=2, lwd=2)
clus <- cutree(hc, k=3)

data <- data.frame(x1 = X[, 1], x2 = X[, 2], clus = factor(clus))

#jpeg('../figures/kmeansResults2.jpeg')
p <- ggplot(data, aes(x1, x2))
p + geom_point(aes(color=clus, shape=clus)) 
#dev.off()


