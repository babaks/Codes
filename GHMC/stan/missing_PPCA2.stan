data {
  int<lower=1> N;  // sample size
  int<lower=1> d;  // dimension of Y_n
  int<lower=1> p;  // lower dimension
  vector[d] Y[N];  // N vector gaussians
  int Nmissing;
  int<lower=1> loc_missing[Nmissing, 2]; // rows and col indicies of missing data 
  int miss[N,d];
  // In R use loc_missing <-which(is.na(Y), arr.ind=TRUE)
}

parameters {
  real<lower=0> sigma2;
  positive_ordered[p] lambda;
  vector[p] u[N];
  vector[d] mu;
  matrix[d,p] X;
}

model {
  lambda ~ cauchy(0,5); 
  sigma2 ~ cauchy(0,5); 
  mu     ~ normal(0,10);
  for (n in 1:N)
    u[n] ~ normal(0,1);
  
  for (r in 1:N) {
    for (c in 1:d) {
      if (miss[r,c]==0)
        Y[r, c] ~ normal(X[c,] * diag_matrix(lambda) * u[r] + mu[c], sigma2);
    }
  }
}

generated quantities{
  vector[Nmissing] predicts;
  for (k in 1:Nmissing) {
      int r;
      int c;
      r <- loc_missing[k, 1];
      c <- loc_missing[k, 2];
      predicts[k] <- X[c,]*diag_matrix(lambda)*u[r] + mu[c];
  }
}

