# This is Simulation 1 in our paper

# The output, which is a 100X4 matrix, is saved as sim_res1.Rdata. Each row shows the area under ROC curve for locFDR, BDP, BRD using v (Eq. 2.7 in the paper) and BRD using average rank respectively. 

library(caTools)
library(fdrtool)

# Number of MCMC iterations and number of burn-in samples
nIter = 2500
burnIn = 500

nData = 100
AUC = matrix(0, nData, 4)

for(d in 1:nData){

# Generating the data
print(d)	
set.seed(d)
z0 = rnorm(460)
z1 = c(rnorm(10, -1), rnorm(10, 1), rnorm(10, 2), rnorm(10, 3))
z = c(z1, z0)

class = c(rep(1, 40), rep(0, 460))


# locFDR
w = fdrtool(z, statistic='normal', plot=FALSE, verbose=FALSE)
fdr = w$qval
p = t(fdr) 
AUC[d, 1] = colAUC(p, class)


# BDP
source('BDP.R')	
res<- cDPODP(z, nIter)
v = NULL
for(i in 1:length(z)){
	
	v[i] = mean(res$matlabel[burnIn:nIter, i] ==1)
	
}
p = t(v) 
AUC[d, 2] = colAUC(p, class)


# BRD
source('BRDZ.R')
simRes = BRD(z, nIter, burnIn)
v1 = NULL
for(i in 1:length(z)){
	
	v1[i] = mean(simRes$post.rank[, i] ==1)
	
}

v2 = colMeans(simRes$post.rank)
	
p = t(v1) 
AUC[d, 3] = colAUC(p, class)

p = t(v2) 
AUC[d, 4] = colAUC(p, class)


save(AUC, file='sim1_res.Rdata')

}