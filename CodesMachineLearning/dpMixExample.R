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
jpeg('../figures/dpData.jpeg')
p <- ggplot(data, aes(x1, x2)) 
p + geom_point()
dev.off()


# Using the example code from the help file in bayesm:

Data1=list(y=X)
Prioralpha=list(Istarmin=1,Istarmax=10,power=.8)
Prior1=list(Prioralpha=Prioralpha)
Mcmc=list(R=1000,keep=1,maxuniq=200)

out1=rDPGibbs(Prior=Prior1,Data=Data1,Mcmc)

rx1=range(X[,1]); rx2=range(X[,2])
x1=seq(from=rx1[1],to=rx1[2],length.out=50)
x2=seq(from=rx2[1],to=rx2[2],length.out=50)
grid=cbind(x1,x2)
plot(out1$nmix,Grid=grid,Data=X)

