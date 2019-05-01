
# GAM
library(gam)
library(MASS)
library(lasso2)

data(Prostate)

# s is for the smoothing spline, we can specify the degrees of freedom as we did for splines. 
gam1 <- gam(lpsa ~ s(age, df=2) + s(lcavol, df=4) + s(lbph, df=3), data = Prostate)

par(ask=T)
plot(gam1, se=TRUE)

summary(gam1)


# Predition; we use the training data for illustration here. In practice, we use the test set as the new dataset. 
y.hat <- predict(gam1, type='response', newdata=Prostate)

plot(Prostate$lpsa, y.hat)


