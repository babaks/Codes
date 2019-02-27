data("metaboliteDataComplete")
mD <- metaboliteDataComplete
Y <- as.matrix(mD)
colnames(Y) <- NULL
rownames(Y) <- NULL
n_data <- dim(Y)[1]
d_obs <- dim(Y)[2]

# compile Stan model for AD
sm  <- stan_model(file.path("stan", "missing_PPCA2.stan"))

d_lat <- 7
perc <- seq(0.1, .8, by=0.1)

for (iii in 1:100) {
for (i in 1:length(perc)) {

  Y_miss <- Y
  for (k in 1:n_data) {
    for (l in 1:d_obs) {
      if (runif(1) < perc[i])
        Y_miss[k,l] <- NA
    }
  }
  loc_miss_binary <- is.na(Y_miss)
  print(dim(loc_miss_binary))
  orig <- Y[loc_miss_binary]
  
  # Run VB PPCA and impute missing
  vb_res <- pca(Y_miss, nPcs=50, method="bpca", verbose=FALSE)
  co <- completeObs(vb_res)
  rec <- co[loc_miss_binary]
  VB_err <- sum(abs(rec-orig)) / sum(loc_miss_binary)
  print(VB_err)
  
  # Grab loadings matrix to find inits for HMC
  L <- loadings(vb_res)
  svdLL    <- svd(L %*% t(L))
  U        <- svdLL$u[, 1:d_lat]
  mu       <- center(vb_res)
  lambda2  <- svdLL$d[1:d_lat]
  
  # sort scales in increasing order
  sorted   <- sort(lambda2, index.return=TRUE)
  U        <- U[,sorted$ix]

  
  Y_miss_trunc <- Y_miss
  Y_miss_trunc[is.na(Y_miss)] = -1
  # HMC Sample
  loc_miss <- which(is.na(Y_miss), arr.ind = TRUE)
  print(dim(is.na(Y_miss)))
  fit <- sampling(sm, data = list(N=n_data,d=d_obs,p=d_lat,Y=Y_miss_trunc, miss=1*is.na(Y_miss),
                                  Nmissing=sum(loc_miss_binary), loc_missing=loc_miss), 
                  iter = 1, chains = 1)
  
  attributes(fit)$inits[[1]]$X <- U 
  attributes(fit)$inits[[1]]$mu <- mu
  
  # Run GHMC
  max_it    <- 1000
  L <- 80
  epsilons=list(sigma2=.002, X=.002, lambda=0.002, mu=.004, u=.004)
  samples <- run_GHMC(fit=fit, max_it=max_it, epsilons=epsilons,
                      L=L, ignore_params=c('u', 'X'), use_stiefel=TRUE,
                      N_generated_quantities = 1, thin=1, burnin=floor(max_it/2))
  
  Nsamples <- length(samples)
  
  #plot(x=1:Nsamples, dist_to_chordal_mean(samples), 'l')
  preds <- laply(samples, .fun=function(x) x$predicts)
  pred <- colMeans(preds)
  HMC_err <- sum(abs(pred-orig)) / sum(loc_miss_binary)
  
  print(HMC_err)
  cat(perc[i],HMC_err, VB_err, "\n",  file='metabolite_errors_per_missing.txt',append=TRUE)
}
}

