# This is Simulation 6 in our paper; for this simulation we use a real dataset: p53. Also, unlike Simulation 1-4, we use the full data as opposed to summary statistics. 

# The output, which is a 100X5 matrix, is saved as sim_res1.Rdata. Each row shows the area under ROC curve for locFDR, BODP, BDP, BRD using v (Eq. 2.7 in the paper) and BRD using average rank respectively. 

library(caTools)
library(fdrtool)

# For this simulation, we use the p53 data set
load('dataP53.Rdata')
Y = data$y
C = unlist(data$c)

n = dim(Y)[1]
nGene = dim(Y)[2]

set.seed(1)
ind = sample(50, 20)
Y = scale(data$y[ind, 1:250])
C = unlist(data$c[ind])
C = sample(C, length(C))

# Number of MCMC iterations and number of burn-in samples
nIter = 2500
burnIn = 500

nData = 100
AUC = matrix(0, nData, 5)

for(d in 1:nData){

# Generating the data
print(d)		
set.seed(d)
y = Y
for(i in 1:5){
	y[C==2, i] = y[C==2, i] + rnorm(1, 0, 1)
}
y = scale(y, scale=FALSE)
	
class = c(rep(1, 5), rep(0, 245))


#locFDR
z = NULL
for(i in 1:(dim(y)[2])){
	z[i] = t.test(y[, i]~C)$p.value	
}

w = fdrtool(z, statistic='pvalue', plot=FALSE, verbose=FALSE)
fdr = w$qval
p = t(fdr) 
AUC[d, 1] = colAUC(p, class)


# BODP
source('BODP.R')
out = BODP(t(y[C==1, ]), t(y[C==2, ]), burnIn, nIter)
odp = out$rd
p = t(odp) 
AUC[d, 2] = colAUC(p, class)


# BDP using full data
source('BDP.R')	
res<- cDPODP_full(y=t(y), tlabel=C, B=nIter, sigma0=1, df0=1, const0 = 2)
v = NULL

for(i in 1:(dim(y)[2])){
	
	v[i] = mean(res$matlabel[burnIn:nIter, i] ==1)
	
}

p = t(v) 
AUC[d, 3] = colAUC(p, class)



# BRD using full data
source('BRDY.R')
simRes <- BRD(y, C, nIter, burnIn)

v1 = NULL
for(i in 1:(dim(y)[2])){
	
	v1[i] = mean(simRes$post.rank[, i] ==1)
	
}

v2 = colMeans(simRes$post.rank)


p = t(v1) 
AUC[d, 4] = colAUC(p, class)

p = t(v2) 
AUC[d, 5] = colAUC(p, class)


print(AUC[d, ])

save(AUC, file='sim6_res.Rdata')

}


