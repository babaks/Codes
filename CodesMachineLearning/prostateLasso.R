library(lars)

#loading prostate data
data <- read.table('prostate.dat', sep=',', header=TRUE)

head(data)

x = (scale(data[, 1:8]))
y = scale(data[, 9], scale=FALSE)

lambda = seq(0.01, 10000, 10)

lm.lasso <- lars(x, y, type='lasso')

plot(lm.lasso)



### We can also use glmnet
library(glmnet)
library(ggplot2)

# Find the optimum lambda using a 10-fold cross validation
cv.lambda <- cv.glmnet(x, y, family='gaussian', nfolds=10)

# Optimum lambda
opt.lambda <- cv.lambda$lambda.min

lasso.res <- glmnet(x, y, family='gaussian', lambda=opt.lambda)

lasso.res$beta

b <- as.vector(lasso.res$beta)
#b[abs(b)<0.02] <- 0
beta <- data.frame(x = seq(1, length(b)), b = b)

ggplot(data=beta, aes(x=x, y=b)) + geom_point(cex=3) + ylab(expression(beta)) + ylim(c(-1.5, 1.5))+
  #theme(plot.margin = unit(c(1, 1, 3, 1),"lines")) + 
  xlab('Variables')+scale_x_discrete(labels=colnames(x)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
