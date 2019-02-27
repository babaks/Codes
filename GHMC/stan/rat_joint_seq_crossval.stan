data {
  int<lower=1> N;  // sample size
  int<lower=1> d;  // dimension of Y_n
  int<lower=1> p;  // dimension of factor matrices
  int Y[N,d];   // N  firings for d neurons

  int Ntrain;
  int trains[Ntrain];  // indicies of training odors
  int tests[N-Ntrain]; // indicies of test odors
  int<lower=0,upper=1> inseq[N];
}

parameters {
  positive_ordered[p] sigma;
  vector[p] u[N];
  vector[d] mu;
  real beta_0;
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
  }
  
  for(n in 1:Ntrain) {
    inseq[trains[n]] ~ bernoulli_logit(beta_0 + dot_product(beta, u[trains[n]]));
  }
}

generated quantities {
  vector[N-Ntrain] inseq_probs;
  vector[N-Ntrain] lpd;
  for (n in 1:(N-Ntrain)) {
    inseq_probs[n] <- inv_logit(beta_0 + dot_product(beta, u[tests[n]]));
    lpd[n] <- bernoulli_logit_log(inseq[tests[n]], beta_0 + dot_product(beta, u[tests[n]]));
  }
}

