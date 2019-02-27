data016.txt: This is the experimental data used in our paper:
Shahbaba, B., Zhou, B., Lan, S., Ombao, H., Moorman, D., and Behseta, S., A Semiparametric Bayesian Model for Detecting Synchrony Among Multiple Neurons, Neural Computation, 26(9), 2025-51.

Main.m: this is the main MATLAB code that calls the sub-functions
CopulaSHMC.m: Performs MCMC calling SphHMC.m, SamplerLatent.m, and SamplerTheta.m
SphHMC.m: Performs spherical HMC sampling
SamplerLatent.m: Performs elliptical slice sampling for latent variables
SamplerTheta.m: Performs slice sampling for hyperparameters
U.m: Calculates the potential and its gradient w.r.t. beta (copula model parameters)
Ut.m: Calculates the potential and its gradient. w.r.t. theta which is the map of beta on standard sphere
LoglikeLatent.m: Calculates the log-likelihood w.r.t. latent variables
LoglikeHyper.m: Calculates the log-likehood w.r.t. hyperparameters
pmfdpmf.m and pmfdpmfDATA.M: Called by u.m to calculate pmf and derivative of pmf
BrownianMat.m: Calculate the covariance matrix of Brownian Motion
BrownianChol.m: Cholesky decomposition of the covariance matrix of Brownian Motion
DetBrownian.m: Calculate the determinant of the covariance matrix of Brownian Motion
SolveBrownian.m: Calculate the inverse of the covariance matrix of Brownian Motion
combrep.m: Derive all possible combinations with replacement


