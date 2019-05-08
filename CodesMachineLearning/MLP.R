rm(list=ls())
library(neuralnet)

set.seed(4)
x = runif(100, -3, 3)
y = (1 + 2*x - 3*sin(x)) + rnorm(100, 0, .5)
plot(x, y)
data <- data.frame(y = y, x = x)
nn <- neuralnet(y~x, data=data, hidden = 5, linear.output = TRUE)

plot(nn)

y.hat <- compute(nn, x)$net.result


pdf('MLP_Example_data.pdf')
plot(x, y)
dev.off()


pdf('MLP_Example_fit.pdf')
plot(x, y)
o <- order(x)
lines(x[o], y.hat[o], lwd=2)
dev.off()





### nnet

library(nnet)
nn <- nnet(x, y, size=50, linout=TRUE, decay = 0.5, data=data)

y.hat <- nn$fitted.values
pdf('MLP_Example_fit_50_decay.pdf')
plot(x, y)
o <- order(x)
lines(x[o], y.hat[o], lwd=2)
dev.off()

