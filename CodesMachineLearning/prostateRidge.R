library(MASS)

#loading prostate data
data <- read.table('prostate.dat', sep=',', header=TRUE)

head(data)

#x = scale(data[data$train==TRUE, 1:8])
#y = scale(data[data$train==TRUE, 9], scale=FALSE)

x =(data[, 1:8])
y =(data[, 9])

x <- as.matrix(x)

# Linear regression model
lm.res <- lm(y~x)
summary(lm.res)


lambda = seq(0.01, 1000, 1)

ridg.res <- lm.ridge(y ~ x, lambda = lambda)

d.f <- NULL
for(i in 1:length(lambda)){

d.f[i] <- sum(diag((x)%*%solve(t(x)%*%(x) + lambda[i]*diag(8))%*%t(x)))	
	
}

beta <- coef(ridg.res)
matplot(d.f, beta[, 2:9], ylim=c(-0.2, 0.8), type='l', xlab=expression(df(lambda)), ylab=expression(beta))

abline(h=0, lty=2, col='gray')

