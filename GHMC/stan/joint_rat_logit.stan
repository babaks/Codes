data {
  int<lower=1> N;  // sample size
  int<lower=1> d;  // dimension of Y_n
  int<lower=1> p;  // dimension of factor matrices
  int Y[N,d];   // N  firings for d neurons

  // multi_logit
  int<lower=0,upper=1> odr[N];
}

parameters {
  positive_ordered[p] sigma;
  real beta_0;
  vector[p] u[N];
  vector[d] mu;
  vector[p] beta; //multi_logit
  matrix[d,p] X;
}

model {
  sigma ~ cauchy(0,5);
  beta_0 ~ normal(0,5);
  beta ~ normal(0,5);
  // u is latent factors
  for (n in 1:N) {
    u[n] ~ normal(0,1);  
  }

  mu ~ normal(0,10); //global mean 
  for(n in 1:N){ //neuron counts and odors
    Y[n] ~ poisson_log(X*(sigma .* u[n]) + mu);
    odr[n] ~ bernoulli_logit(dot_product(beta, u[n])+beta_0);
  }
}

