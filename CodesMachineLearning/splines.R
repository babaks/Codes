library(splines)

set.seed(4)
x = runif(100, -3, 3)
y = (1 + 2*x - 3*sin(x)) + rnorm(100, 0, .5)
plot(x, y)

abline(v=-1, col='black')
abline(v=1, col='black')


# Piecewise constant
h1 = x <=- 1
h2 = x > -1 & x <= 1
h3 = x > 1 

H = cbind(h1, h2, h3) 

beta = solve(t(H)%*%H)%*%t(H)%*%y

y.hat = H%*%beta

lines(c(-3, -1), c(beta[1], beta[1]), lwd=2)
lines(c(-1, 1), c(beta[2], beta[2]), lwd=2)
lines(c(1, 3), c(beta[3], beta[3]), lwd=2)


# Piecewise linear
h1 = x <= -1
h2 = x > -1 & x <= 1
h3 = x > 1 
h4 = h1*x
h5 = h2*x
h6 = h3*x

H = cbind(h1, h2, h3, h4, h5, h6) 

beta = solve(t(H)%*%H)%*%t(H)%*%y

y.hat = H%*%beta

plot(x, y)
abline(v=-1, col='black')
abline(v=1, col='black')

lines(c(-3, -1), c(beta[1]+beta[4]*(-3), beta[1]+beta[4]*(-1)), lwd=2)
lines(c(-1, 1), c(beta[2]+beta[5]*(-1), beta[2]+beta[5]*(1)), lwd=2)
lines(c(1, 3), c(beta[3]+beta[6]*(1), beta[3]+beta[6]*(3)), lwd=2)




# Piecewise linear with continuity
h1 = 1
h2 = x
h3 = pmax(0, x - (-1))
h4 = pmax(0, x - (1))

H = cbind(h1, h2, h3, h4) 

beta = solve(t(H)%*%H)%*%t(H)%*%y

y.hat = H%*%beta

plot(x, y)
abline(v=-1, col='black')
abline(v=1, col='black')

lines(c(-3, -1), c(beta[1]+beta[2]*(-3), beta[1]+beta[2]*(-1)), lwd=2)
lines(c(-1, 1), c(beta[1] - beta[2], beta[1]+beta[2] + 2*beta[3]), lwd=2)
lines(c(1, 3), c(beta[1]+beta[2] + 2*beta[3], beta[1]+3*beta[2]+4*beta[3]+2*beta[4]), lwd=2)



# If we use B-splines

plot(x, y)

abline(v=-1, col='black')
abline(v=1, col='black')

basis.fn <- bs(x, knots=c(-1, 1), degree=1)
B.model <- lm(y~basis.fn)
beta = coef(B.model)

t <- seq(-3, 3, .1)
basis.fn.t = predict(basis.fn, t)
y.hat <-cbind(1, basis.fn.t)%*%beta
lines(t, y.hat, lwd=2)




# Order-4 polynomial spline (degree=3)

plot(x, y, xlim=c(-5, 5), ylim=c(-10, 10))

abline(v=-1, col='black')
abline(v=1, col='black')

basis.fn <- bs(x, knots=c(-1, 1), degree=3)
cubicSp.model <- lm(y~basis.fn)
beta = coef(cubicSp.model)

t <- seq(-5, 5, .1)
basis.fn.t = predict(basis.fn, t)
y.hat <-cbind(1, basis.fn.t)%*%beta
lines(t, y.hat, lwd=2)


# Natural cubic spline

basis.fn <- ns(x, knots=c(-1, 1))
beta = coef(lm(y~basis.fn))

t <- seq(-5, 5, .1)
basis.fn.t = predict(basis.fn, t)
y.hat <-cbind(1, basis.fn.t)%*%beta
lines(t, y.hat, lwd=2, lty=2)


# Smoothing splines
plot(x, y)

abline(v=-1, col='black')
abline(v=1, col='black')

ss.model <- smooth.spline(x, y, spar=0.1)
lines(ss.model, lwd=2)

ss.model <- smooth.spline(x, y, spar=0.9)
lines(ss.model, lwd=2)


ss.model <- smooth.spline(x, y, df = 30)
lines(ss.model, lwd=2)

ss.model <- smooth.spline(x, y, df = 10)
lines(ss.model, lty=2, lwd=2)

ss.model <- smooth.spline(x, y, df = 3)
lines(ss.model, lty=3, lwd=2)


legend('topleft', legend=c(expression(paste(df[lambda], '=30')), expression(paste(df[lambda], '=10')), expression(paste(df[lambda], '=3')) ), lty=c(1, 2, 3), lwd=2)

