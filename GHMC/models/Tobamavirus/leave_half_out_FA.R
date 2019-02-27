load("data/tobamovirus.Rdata")
Y <- as.matrix(virusdat)
colnames(Y) <- NULL

#compile
sm  <- stan_model('stan/Grassmann_FA.stan')

cond_mean <- function(obs, mu, Sigma, indices){
  S_12 <- Sigma[indices, -indices]
  S_22 <- Sigma[-indices, -indices]
  cm <- mu[indices] + as.vector(S_12 %*% solve(S_22, (obs[-indices]- mu[-indices])))
  return(cm)
}


######################################################################
#start loho
Nreps <- 100
p <- 2
N <- 38
MLE_error <- array(0,dim=c( N,16,Nreps))
BAY_error <- array(0,dim=c( N,16,Nreps))

par(mfrow=c(2,2), mar=c(2,2,2,2))
for(j in 1:N) {

  Y_j     <- Y[j,]
  Y_minus <- Y[-j,]
  S        <- Y_minus-matrix(rep(colMeans(Y_minus), N-1), nrow=N-1, byrow=TRUE)
  fa       <- factanal(S, factors=p, rotation="none")
  sigma2   <- fa$uniquenesses
  tau_est  <- sqrt(aaply(S, .fun=var, .margins=2))
  Sig      <- fa$loadings %*% t(fa$loadings)
  U        <- svd(Sig)$u[,1:p]
  C        <- diag(tau_est) %*% Sig %*% diag(tau_est) + diag(sigma2)

  # run HMC
  fit <- sampling(sm, data = list(N=N-1,d=18,p=p,Y=Y_minus), 
                  iter = 1, chains = 1)
  
  attributes(fit)$inits[[1]]$tau <- tau_est
  attributes(fit)$inits[[1]]$sigma2 <- sigma2 
  attributes(fit)$inits[[1]]$X <- U 
  attributes(fit)$inits[[1]]$mu <- colMeans(Y_minus) 
  max.it    <- 200
  L       <- 50
  epsilons <- list(X=.02, sigma2=0.01, tau=0.01, mu=0.02)
  samples <- run_GHMC(fit=fit, max_it=max.it, 
                      epsilons=epsilons,
                      L=L, ignore_params=c(), use_stiefel=FALSE,
                      N_generated_quantities = 1, thin=1, burnin=0)
  Nsamples <- length(samples)
  
  plot(x=1:Nsamples, dist_to_chordal_mean(samples, X0=U), 'l')
  for(k in 1:16){
    for(m in 1:Nreps){
      indices <- sample(1:18,k)
      #get conditional means from MLE
      Yhat <- cond_mean(obs=Y_j, mu=colMeans(Y_minus),
                            Sigma=C, indices=indices)
      MLE_error[j,k,m] <- sum(abs(Yhat-Y_j[indices]))/k
    
      #get conditional mean for ith iteration of MCMC (from 1001 to max)
      Yhat_bay <- matrix(0,k,Nsamples/2)
      for(i in (Nsamples/2+1):Nsamples){
        s <- samples[[i]]
        sigma_cov <- s$C
        Yhat_bay[,i-max.it/2] <- cond_mean(obs=Y_j, mu=s$mu,
                              Sigma=sigma_cov, indices=indices)
       }
      BAY_error[j,k,m] <- sum(abs(rowMeans(Yhat_bay)-Y_j[indices]))/k
    }
  }
  print(cbind(rowMeans(BAY_error[j,,], dims=1), rowMeans(MLE_error[j,,], dims=1)))
  
}#end LOO

#print errors

save(BAY_error, MLE_error, file="tobamo_FA_errs_grass.Rdata")









