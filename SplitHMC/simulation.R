# This file runs the simulation example by calling logisticHMC.R

example = 'simulatedData'

# Generate Data
set.seed(2011)
beta = rnorm(101, 0, 1)
n = 10000
Sigma = matrix(0, 100, 100)
diag(Sigma) <- c(rep(25, 5), rep(1, 5), rep(0.2^2, 90))
X = rmvnorm(n, sigma=Sigma)
x =  cbind(1, X)
p = exp(x%*%beta)/(1+exp(x%*%beta))
y = rbinom(n, 1, prob=p)


# Set the MCMC parameters
leapFrog = 20
epsilon = 0.015

leapFrog.splitNorm = 10
epsilon.splitNorm = (leapFrog*epsilon)/leapFrog.splitNorm

leapFrog.splitData = 3
epsilon.splitData = (leapFrog*epsilon)/leapFrog.splitData
M = 9
frac = 0.6


# Number of Iterations
nIter = 50000

# Standard deviation of regression parameters. In practice, we usually use hyperpriors for sigma.
sigma = 5

# Create a folder called "Results" to save the final results
dir.create('./Results', showWarnings = FALSE)


# Obtain all the required functions for running the codes. Also, logisticHMC.R finds MAP and Fisher Information. Additionally, it splits the data to (x0, y0) for U0 and (x1, y1) for U1
source('logisticHMC.R')


# The following codes run the models and measure the performances of the three methods (Standard HMC, Split HMC with normal approximation, Split HMC with data splitting) in terms of the CPU time per iteration, acceptance probability (acp), autocorrelation time (act), and efficiency, which is defined as the product of CPU and act and is interpreted as the amount of time required to genereate an independent sample. The program creates a folder called "Results" in the current directory and saves the results in three .Rdata files (one for each method) there.  

# Run Standard HMC
start.time = proc.time()
samp <- hmcStandard(x, y, q.hat, nIter, epsilon, leapFrog)
end.time = proc.time()
cpu.time = end.time[1] - start.time[1]
results <- getPerformance(samp, cpu.time, nIter)
save(samp, results, epsilon, leapFrog, file=paste('Results/', example, '_hmcStandard.Rdata', sep = ''))


# Run Split HMC with normal approximation
start.time = proc.time()
samp <- hmcSplitNorm(x, y, q.hat, FI, R, nIter, epsilon.splitNorm, leapFrog.splitNorm)
end.time = proc.time()
cpu.time = end.time[1] - start.time[1]
results <- getPerformance(samp, cpu.time, nIter)
save(samp, results, epsilon.splitNorm, leapFrog.splitNorm, file=paste('Results/', example, '_hmcSplitNorm.Rdata', sep = ''))


# Run Split HMC with data splitting
start.time = proc.time()
samp <- hmcSplitData(x, y, x0, y0, x1, y1, q.hat, nIter, epsilon.splitData, leapFrog.splitData, M)
end.time = proc.time()
cpu.time = end.time[1] - start.time[1]
results <- getPerformance(samp, cpu.time, nIter)
save(samp, results, epsilon.splitData, leapFrog.splitData, file=paste('Results/', example, '_hmcSplitData.Rdata', sep = ''))
	
