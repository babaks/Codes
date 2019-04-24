library(MASS)
library(lasso2)
library(pls)
library(glmnet)

data(Prostate)

head(Prostate)

Prostate$psa <- exp(Prostate$lpsa)

# attach(Prostate)


# Hypothesis testing; H: seminal vesicle invasion is associated with prostate specific antigen after controlling for the age of patients. 

lm1 <- lm(psa ~ svi + age, data = Prostate)
summary(lm1)

# How does the interpretation of the model parameters change if we use lpsa

lm2 <- lm(lpsa ~ svi + age, data = Prostate)
summary(lm2)


# Prediction; we want to predict psa given the predictors. 

Prostate$psa <- NULL

lm3 <- lm(lpsa ~ ., data = Prostate)
summary(lm3)

# What is wrong with using the same data to evaluate our model?
lpsa.hat <- predict(lm3, new.data=Prostate)
mean((Prostate$lpsa - lpsa.hat)^2)

# Data splitting
n = dim(Prostate)[1]
id.train <- sample(n, round(n*0.6))
id.remain <- setdiff(1:n, id.train)
id.val <- sample(id.remain, round(n*0.2))
id.test <- setdiff(id.remain, id.val)

# If we don't have fine-tuning of parameters, we can combine training and validation
lm4 <- lm(lpsa ~ ., data = Prostate[c(id.train, id.val), ])
x.new <- as.data.frame(Prostate[id.test, 1:8])
lpsa.hat <- as.matrix(cbind(1, x.new)) %*% matrix(coef(lm4))
mean((Prostate$lpsa[id.test] - lpsa.hat)^2)


# Instead of using all the variables, we might want to use the best subset, for example, based on AIC. 
# NOTE that as discussed in the class, this is just an exploratory method if we don't have a well-defined hypothesis. 

lm5 <- stepAIC(lm3, scope = list(lower = ~1), direction='both')
summary(lm5)

AIC(lm3)
AIC(lm5)


# Using BIC
lm6 <- stepAIC(lm3, scope = list(lower = ~1), direction='both', k = log(nobs(lm3)))
summary(lm6)

BIC(lm3)
BIC(lm6)




Prostate[, 1:8] <- scale(Prostate[, 1:8])

x <- as.matrix(Prostate[, 1:8])
y <- Prostate[, 9]
n <- length(y)

# PCR

s <- t(x)%*%x/n
e <- eigen(s)
lam <- e$values
v <- e$vectors

plot(lam, type='l') # scree plot

z <- x%*%v

# The derived variable are ordered according to their sample variance
apply(z, 2, var)

# We could also use princomp

pc <- princomp(x)
z<- pc$scores[, 1:5] # Choosing the first 5 principal components 
pcr.mod <- lm(y ~ z)


# We could have the data-splitting method to choose the number of principal components 



# Here, we use 5-fold CV to evaluate the performance of the PCR with 5 components

# To be even more precise, we should put this within the loop too
pc <- princomp(x)
z<- pc$scores[, 1:5] # Choosing the first 5 principal components 

ind <- sample(5, n, replace=TRUE)


mse = rep(0, 5) 
for(i in 1:5){

	pcr.mod <- lm(y[ind!=i] ~ z[ind!=i, ])

	y.hat <- as.matrix(cbind(1, z[ind==i, ])) %*% matrix(coef(pcr.mod))
	
	mse[i] = mean((y[ind==i] - y.hat)^2)
	
	
} 

mean(mse)



# we can use the pcr function in the pls package

pcr.mod <- pcr(y ~ x, ncomp=5, scale=TRUE, validation = "CV", segments = 5)


# PLS

pls.mod <- plsr(y ~ x, ncomp=5, data=Prostate, scale=TRUE)


# We can use the "validation" for fine-tuning parameters. See the help document for plsr.

# Ridge regression

ridge.res <- lm.ridge(y~x, lambda = seq(0, 10, 0.1))
plot(ridge.res)
select(ridge.res)

# Lasso with glmnet

lasso.fit=glmnet(x,y)
plot(lasso.fit)

y.hat <- predict(lasso.fit, x, s=0.01, type='response')


