# This file runs the chess example by calling logisticHMC.R

example = 'chess'

# Get the data
data = read.table('chess.txt', header = FALSE, sep=',')
y = as.numeric(data[, 37]=='won')
X = data[, 1:36]
x = matrix(0, dim(X)[1], 36)
for(i in 1:36){
	x[, i] = as.numeric(X[, i])-1
}
x =  cbind(1, x)

# Set the MCMC parameters
leapFrog = 20
epsilon = 0.09

leapFrog.splitNorm = 9
epsilon.splitNorm = (leapFrog*epsilon)/leapFrog.splitNorm

leapFrog.splitData = 2
epsilon.splitData = (leapFrog*epsilon)/leapFrog.splitData
M = 15
frac = 0.65

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
#results <- getPerformance(samp, cpu.time, nIter)
#save(samp, results, epsilon, leapFrog, file=paste('Results/', example, '_hmcStandard.Rdata', sep = ''))


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
	

