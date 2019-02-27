data {
  int<lower=1> N;  // sample size
  int<lower=1> d;  // dimension of Y_n
  int<lower=1> p;  // lower dimension
  vector[d] Y[N];   //N vector gaussians
}

parameters {
  vector<lower=0>[d] sigma2;
  vector<lower=0>[d] tau;
  vector[d] mu;
  matrix[d,p] X;
}

model {
  tau ~ cauchy(0,5); 
  sigma2 ~ cauchy(0,5); 
  mu     ~ normal(0,10);
  {
    matrix[d,d] L;
    L <- cholesky_decompose(quad_form_diag(X * (X') + diag_matrix(sigma2), tau));
  
    for(n in 1:N){
      Y[n] ~ multi_normal_cholesky(mu, L); 
    }
  }
}

generated quantities {
  matrix[d,d] C;
  C <- quad_form_diag(X * (X') + diag_matrix(sigma2), tau);
}
