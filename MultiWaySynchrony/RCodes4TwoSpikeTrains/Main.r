
#######################################################################################################################
## Perform analysis of spike train data
#######################################################################################################################

load("data016.RData") ## Load Data

## Source functions to do MCMC sampling
source( "sampler.latent.r" )
source( "sampler.hyper.r" )
source( "MCMC.r" )
source( "log.like.r" )
source( "solve.brownian.r" )
source( "BrownianMat.r" )

## Source packages required
library( splines )
library( spatstat )
library( coda )
library(compositions)


## Initialization for MCMC sampling
set.seed(123456)
time <- 5*(0:(dim(N11)[1]-1))/1000 
Niteration <- 3000
BurnIn <- 1000
w <- 1
m <- 10
epsilon <- 0.02
Bound <- 20

## Analyze pairs of neurons( Neuron 1 and Neuron 7) under lever 1
x <- time[1:2000]
t1 <- N1[1,1:2000,]
t2 <- N1[7,1:2000,]
Sample <- BM_LAG_TTRE( x, t1, t2, Niteration, w, m, epsilon, Bound )
save.image( "data016_N1N7_lever1.RData" )

## Analyze pairs of neurons( Neuron 1 and Neuron 7) under lever 2
x <- time[1:2000]
t1 <- N2[1,1:2000,]
t2 <- N2[7,1:2000,]
Sample <- BM_LAG_TTRE( x, t1, t2, Niteration, w, m, epsilon, Bound )
save.image( "data016_N1N7_lever2.RData" )

## Analyze pairs of neurons( Neuron 7 and Neuron 9) under lever 1
x <- time[1:2000]
t1 <- N1[7,1:2000,]
t2 <- N1[9,1:2000,]
Sample <- BM_LAG_TTRE( x, t1, t2, Niteration, w, m, epsilon, Bound )
save.image( "data016_N7N9_lever1.RData" )

## Analyze pairs of neurons( Neuron 7 and Neuron 9) under lever 2
x <- time[1:2000]
t1 <- N2[7,1:2000,]
t2 <- N2[9,1:2000,]
Sample <- BM_LAG_TTRE( x, t1, t2, Niteration, w, m, epsilon, Bound )
save.image( "data016_N7N9_lever2.RData" )



