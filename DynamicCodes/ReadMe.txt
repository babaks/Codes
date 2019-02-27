Data016_100ms.txt: This is the experimental data used in our paper.

Main.m: This is the main function that calls all the other functions required for inference. 

MCMC.m: This function Performs MCMC by calling subfunctions to sample the posterior distribution of latent variables, thresholds, hyperparameters and so forth.

samplerLatent.m: This function draws samples for latent variables from truncated normal distributions.

samplerTr.m: This function draws samples for thresholds using elliptical slice sampler.

samplerEta.m, samplerRho.m, samplerAlpha.m and samplerJ2.m: These four functions are used to sample the hyperparameters of Gaussian Process, respectively.

samplerSigma.m: This function draws samples for covariance matrix by using ellptical slice sampler.

samplerInd.m: This function draws samples for the indicator variables.

samplerQt.m: This function draws samples for q_t which is the latent process for the indicators.

samplerTheta.m: This function draws samples for the hyperparameters of the latent process, q_t.

laprnd.m: This function draws samples from laplace distribution.

rmvnrnd.m: This function draws samples from multivariate truncated normal distributions.

TruncatedGaussian.m: This function draws samples from univariate truncated normal distributions.

loglikeTr.m: This function calculates the log-likelihood function w.r.t. the threshods.
loglikeBmHyper.m: This function calculates the log-likelihood function w.r.t. the hyperameters of brownian motion.
loglikeGpHyper.m: This function calculates the log-likelihood function w.r.t. the hyperameters of Gaussian Process.
loglikeQt.m: This function calculates the log-likelihood function w.r.t. the latent process, q_t, of the indicators.
loglikeSigma.m: This function calculates the log-likelihood function w.r.t. the covariance matrix.



