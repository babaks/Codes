# This is Simulation 4 in our paper; for this simulation we use a real dataset: p53.

# The output, which is a 100X4 matrix, is saved as sim_res1.Rdata. Each row shows the area under ROC curve for locFDR, BDP, BRD using v (Eq. 2.7 in the paper) and BRD using average rank respectively. 

library(caTools)
library(fdrtool)

# For this simulation, we use the p53 data set
load('dataP53.Rdata')
Y = data$y
C = unlist(data$c)

nData = 100
n = dim(Y)[1]
nGene = dim(Y)[2]

Z = NULL
for(i in 1:nGene){

	Z[i] = t.test(Y[, i]~C)$statistic 
	
}

ord = order(abs(Z), decreasing=TRUE)
Z = Z[ord]




# Number of MCMC iterations and number of burn-in samples
nIter = 2500
burnIn = 500

nData = 100
AUC = matrix(0, nData, 4)

for(d in 1:nData){

# Generating the data
print(d)	
set.seed(d)
z0 = sample(Z[501:nGene], 460)
z1 = c(sample(Z[1:100], round(runif(1, 5, 10))), sample(Z[201:400], round(runif(1, 20, 30))))
z = c(z1, z0)

class = c(rep(1, length(z1)), rep(0, length(z0)))


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


save(AUC, file='sim4_res.Rdata')

}