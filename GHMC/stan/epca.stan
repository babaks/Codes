data {
  int<lower=1> N;  // sample size
  int<lower=1> d;  // dimension of Y_n
  int<lower=1> p;  // dimension of factor matrices
  int<lower=-1, upper=1> Y[N,d];   // N  firings for d neurons
  int<lower=0, upper=1> Ytrue[N,d];   // N  firings for d neurons
  int miss[N,d]; // which entries are missing
}

parameters {
  positive_ordered[p] lambda;
  vector[p] u[N];
  vector[d] mu;
  matrix[d,p] X;
}

model {
  lambda ~ cauchy(0,.5);
  for (n in 1:N) {
    // latent factors
    u[n] ~ normal(0,lambda);  
  }

  mu ~ normal(0,1); //global mean 
  for(n in 1:N){ 
    for (m in 1:d) {
      if (miss[n,m] == 0)
        Y[n,m] ~ bernoulli_logit(X[m,]*(u[n]) + mu[m]);
    }
  }
}

generated quantities {
  real log_lik_miss;
  matrix[N,d] predicts;
  real alllik;
  alllik <- cauchy_log(lambda, 0, .5) + normal_log(mu, 0, 1);
  
  log_lik_miss <- 0;
  for (n in 1:N) {
    alllik <- alllik + normal_log(u[n], 0, lambda); // + bernoulli_logit_log(Y[n,], X*u[n] + mu);
    for (j in 1:d) {
      real ll;
      ll <- bernoulli_logit_log(Ytrue[n,j], X[j,]*(u[n]) + mu[j]);
      if (miss[n,j] == 1)
        log_lik_miss <- log_lik_miss + ll;
      alllik <- alllik + ll;
      predicts[n, j] <- bernoulli_logit_log(1, X[j,]*(u[n]) + mu[j]);
    }
  }
}

