# We assume the price has a normal(mu, sigma^2) distribution
mu=30
sigma = 5

x = seq(mu-3*sigma, mu+3*sigma, .01)
expected.gain = x*dnorm(x, mu, sigma)
plot(x, y)

# Setting the derivative to zero
print(0.5*(mu+sqrt(mu^2+4*sigma^2)))

# Finding the maximum value of expected gain numerically
answer = x[which.max(y)]
print(answer)

abline(v=mu, lty=1)
abline(v=answer, lty=2)
