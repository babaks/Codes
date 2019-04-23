x <- USArrests
x <- scale(x)

#PCA
pca.res = princomp(x)
summary(pca.res)
loadings(pca.res) # small values are not printed
z <- pca.res$scores # scores 

# alternative view based on SVD:
s <- svd(x)
s$v

# Note that we can recover x as follows:
x.hat <- z%*%t(s$v)


pdf('pca_usarrests.pdf')
biplot(pca.res)
dev.off()


# since the above dataset is low dimensional, we continue our illustration using the state.x77 dataset available in R
x <- scale(state.x77)
pca.res <- princomp(x)

pdf('pca_states.pdf')
biplot(pca.res)
dev.off()

#scree plot
plot(pca.res)
df <- data.frame(PC = 1:length((pca.res$sdev)), Variance = (pca.res$sdev)^2)
pdf('scree_states.pdf')
qplot(PC, Variance, data= df, geom='line')
dev.off()
# cumulative percent of variance explained

cumvar <- 100*stats:::summary.prcomp(pca.res)$importance[3, ]
df <- data.frame(PC = 1:length((cumvar)), VarianceExplained = cumvar)

pdf('var_states.pdf')
qplot(PC, cumvar, data=df, geom='line', ylab='Cumulative variance (%)')
dev.off()

#Scores
z = pca.res$scores[, 1:4]



#Factor Analysis
fa.res <- factanal(x, 2, scores='regression')
print(fa.res)

pdf('fa_states.pdf')
biplot(fa.res$scores[, 1:2], fa.res$loadings[,1:2] )
dev.off()


# ICA
library(fastICA)

# Examples from the fastICA package

# Example1: data generated from two independent uniform variables
S <- matrix(runif(1000), 500, 2)

# A is the mixing matrix
A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)

# This is the observed data
X <- S %*% A

# using ICA to identify two independent components
a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
  method = "C", row.norm = FALSE, maxit = 200,
  tol = 0.0001, verbose = TRUE)

pdf('ica_example1.pdf')
par(mfrow = c(2, 2))
plot(S, main = "Original data",
  xlab = "", ylab = "")
plot(a$X, xlab = "", ylab = "", main = "Pre-processed data")
plot(a$X %*% a$K,xlab = "", ylab = "", main = "PCA components")
plot(a$S, xlab = "", ylab = "", main = "ICA components")
dev.off()


# Example 2: two independent sources
S <- cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A
a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
  method = "R", row.norm = FALSE, maxit = 200,
  tol = 0.0001, verbose = TRUE)

pdf('ica_example2.pdf')
par(mfcol = c(2, 4))
plot(1:1000, S[,1 ], type = "l", main = "Original Signals",
  xlab = "", ylab = "")
plot(1:1000, S[,2 ], type = "l", xlab = "", ylab = "")
plot(1:1000, X[,1 ], type = "l", main = "Mixed Signals",
  xlab = "", ylab = "")
plot(1:1000, X[,2 ], type = "l", xlab = "", ylab = "")
plot(1:1000, (a$X %*% a$K)[, 1], type = "l", main = "PCA source estimates",
  xlab = "", ylab = "")
plot(1:1000, (a$X %*% a$K)[, 2], type = "l", xlab = "", ylab = "")
plot(1:1000, a$S[,1 ], type = "l", main = "ICA source estimates",
  xlab = "", ylab = "")
plot(1:1000, a$S[, 2], type = "l", xlab = "", ylab = "")
dev.off()

