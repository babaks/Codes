data016.RData: This is the experimental data used in our paper.

Main.m: This is the main R code that calls the sub-functions to analyze spike train data
MCMC.R: Performs MCMC calling sampler.hyper.R and sampler.latent.R
sampler.latent.R: Sample latent variables using elliptical slice sampling, sample extra terms using Metropolis Hastings, 
sample lag from multinomial distribution and sample prior probabilities of lags from Dirichlet distribution.
sampler.hyper.R: Sample hyperparameters using slice sampler.
solve.brownian.R: Calculate the inverse and the determinant of the covariance matrix of Brownian Motion 
BrownianMat.m: Calculate the covariance matrix of Brownian Motion and its Cholesky decomposition.

