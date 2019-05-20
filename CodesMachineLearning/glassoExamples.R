library(glasso)
library(glmnet)
library(network)
library(GGally)
library(sna)
library(ggplot2)

# continuous data
data <- read.csv('NeuronsContinuous.csv')

s<- var(data)

# Graphical lasso
a<-glasso(s, rho=2)

# sparse inverse covariance
P <- a$wi
A <- ifelse(P!=0 & row(P)!=col(P),1,0)
g <- network(A, directed = FALSE)
network.vertex.names(g) = paste('Neuron', 1:6, sep='')

jpeg('Neuron_MorkovGraph_Continuous.jpeg')
ggnet2(g, size = 18, label = TRUE, label.size = 2.5)
dev.off()


# Discrete data
data <- read.csv('NeuronsDiscrete.csv')

# Logistic regression coefficients
beta <- matrix(0, 6, 6)
for (j in 1:6){
  
  y <- data[, j]
  X <- do.call(cbind, data[, -j])
  
  mod <- glmnet(X, y, family='binomial', lambda=0.1)
  
  b <- as.vector(mod$beta)
  beta[j, -j] <- b
  
}

# min symmetrization
for(j in 1:6){
  for(k in j:6){
    beta[j, k] <- min(abs(beta[j, k]), abs(beta[k, j]))
    beta[k, j] <- beta[j, k]
  }
}

g <- network(beta, directed = FALSE)
network.vertex.names(g) = paste('Neuron', 1:6, sep='')

pdf('Neuron_MorkovGraph_Discrete.pdf')
par(mar=c(3, 3, 3, 3)) 
ggnet2(g, size = 18, label = TRUE, label.size = 2.5)
dev.off()

