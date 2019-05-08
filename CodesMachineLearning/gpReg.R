library(mvtnorm)


# Sampling from the prior on y(x)

lambda = 2
eta = 1
rho = 0.5
sigma = 0.1

n = 200
x = runif(n, -5, 5)

diffMatAll = matrix(x, nrow=n, ncol=n) - matrix(x, nrow=n, ncol=n, byrow=TRUE)

# Covariance matrix

# Linear regression model
# C = eta + rho*(x%*%t(x))+ sigma*diag(1, nrow=n, ncol=n)

# Nonlinear regression model
C = lambda*exp(-(diffMatAll^2)) + eta*exp(-rho*(diffMatAll^2))+ sigma *diag(1, nrow=n, ncol=n)

# Nonlinear regression model with a different covariance function
# C = lambda + eta*exp(-rho*(diffMatAll^2))+ sigma*diag(1, nrow=n, ncol=n)

y = rmvnorm(1, sigma=C)

plot(x, y)

# Generating data

n = 200
set.seed(2)
x = runif(n, -5, 5)
y = (1 -.6*x + .1*x^2 - 3.5*sin(x) + 2.2*cos(x)) + rnorm(n, 0, 1)
#y = .5*(x^3) - 3*x + 1 + rnorm(n, 0, .5)
plot(x, y, xlim=c(-5, 5))


# Dividing the data into training and test

ind.tr <- sample(n, 50)
ind.te <- setdiff(seq(1, n), ind.tr)

x.tr <- x[ind.tr]
y.tr <- y[ind.tr]
x.te <- x[ind.te]
y.te <- y[ind.te]

nTrain = length(x.tr);
nTest = length(x.te);

# Function to get the posterior predictive probabiltiy
gpReg = function(x.tr, y.tr, x.te, eta, rho, sigma){

x = c(x.tr, x.te);

n = length(x);


diffMatAll = matrix(x, nrow=n, ncol=n) - matrix(x, nrow=n, ncol=n, byrow=TRUE)

C = lambda + eta*exp(-rho*(diffMatAll^2))+ sigma*diag(1, nrow=n, ncol=n);

Ctrn = C[1:nTrain, 1:nTrain];
invCtrn = solve(Ctrn)

K = C[1:nTrain, (nTrain+1):n];
v = C[(nTrain+1):n, (nTrain+1):n];

# E(y.te | y.tr)
y.hat = t(K)%*%invCtrn%*%y.tr;

# Var(y.te | y.te)
v.hat = v - t(K)%*%invCtrn%*%K; 

return(list(y.hat=y.hat, v.hat=v.hat))

}


res = gpReg(x.tr, y.tr, x.te, eta, rho, sigma)

y.hat = res$y.hat
v.hat = res$v.hat

mse = mean((y.te - y.hat)^2)



# plotting the model over a grid
x.te <- seq(-5, 5, .1)
res = gpReg(x.tr, y.tr, x.te, eta, rho, sigma)
y.hat = res$y.hat
v.hat = res$v.hat

lines(x.te, y.hat, lwd=2)



### You can also use the function gausspr in kernlab

library(kernlab)

data.tr <- data.frame(x = x.tr, y = y.tr)
data.te <- data.frame(x = x.te, y = NA)

gpm <- gausspr(y~x, data=data.tr)
y.hat <- predict(gpm, newdata=data.te)
lines(x.te, y.hat, col='red', lwd=2)



# We can also use this function for classification 
data = read.table('Pima.txt', header=TRUE)
n = dim(data)[1]
p = dim(data)[2]

ind.tr <- sample(n, floor(2*n/3))
ind.te <- setdiff(seq(1, n), ind.tr)

gpm <- gausspr(diabetes~ ., type = 'classification', kernel='rbfdot', data=data[ind.tr, ])
pred.class <- predict(gpm, data[ind.te, ])
act.class <- data$diabetes[ind.te]
table(act.class, pred.class)
mean(act.class == pred.class)
