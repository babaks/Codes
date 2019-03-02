The files in this folder are MATLAB and R programs based on the examples discussed in the STATS 225 course. This is the list of these programs:

rejectionSampling: This program performs rejection sampling from Beta(3, 10) using a Uniform(0,1).

goodImportanceSampling: This is an example of a good importance sampling. Here we want to use t_3 for sampling from N(0, 1).

badImportanceSampling: Here, we are using N(0, 1) to sample from t_3. This is an example where the importance sampling will not work properly. 

continuousMC1: This is a continuous space Markov chain with N(x/2, sqrt(3/4)) transition distribution.

binoPostMetropolis: This program uses the Metropolis algorithm for the binomial model with Beta prior.

poissonPostMH: This program uses the Metropolis-Hastings (MH) algorithm for the Poisson model with Gamma prior. 

normPostMetropolis1: This program uses the Metropolis algorithm for the multivariate normal model with unknown mean and known covariance matrix. All parameters are updated simultaneously. 

normPostMetropolis2:  This program uses the Metropolis algorithm for the multivariate normal model with unknown mean and known covariance matrix. The parameters are updated one at a time. 

normPostGibbs: This program uses the Gibbs sampler for simulating from the posterior distribution of a normal model with unknown mean and variance

modelHousePrice: This program uses the Gibbs sampler for the hierarchical model of US house pricing described in the session 9.

modelRats: This program models rats tumor using a hierarchical Bayesian model.

ratTumor.txt: This is the dataset for the rats tumor model. 

BayesianLogitSnore: This program uses logistic regression to model the effect of snoring on heat disease.

linRegKidsScore: This is a linear regression model applied to the kids score dataset.

linRegKidsScoreNonInfo: This is a linear regression model with noninformative priors applied to the kids score dataset.

kidiq.txt: This is the kids score dataset used in the linRegKidsScore and linRegKidsScoreNonInfo.

BayesianLogitTitanic: This program uses logistic regression to model the effect of social class in Titanic survival.

titanic.dat: This is the titanic dataset obtained from the delve website (http://www.cs.utoronto.ca/~delve/). 

BayesianLogitBCW: This program uses logistic regression model to predict the risk of breasst cancer.  

bcw.data: This the Wisconsin breast cancer data obtained from the UCI repository (http://mlearn.ics.uci.edu/MLRepository.html).

