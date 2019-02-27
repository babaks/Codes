library(caret)

######################################################################
SuperChris1 <- read.csv(file.path("data", "counts_superchris_1.csv"), header=FALSE)
Y <- as.matrix(SuperChris1)
Y <- Y

inseq <- read.csv(file.path("data", "seq_superchris_1.csv"), header=FALSE)
inseq   <- as.matrix(inseq)
colnames(inseq) <- c()
inseq   <- as.vector(inseq)
N     <- dim(Y)[1]
d_obs <- dim(Y)[2]
d_lat <- 2

# compile Stan model for AD
sm  <- stan_model(file.path("stan", "rat_joint_seq_crossval.stan"))

# Cross validation splitting into Kfolds
Nfolds <- 10
flds <- createFolds(inseq, k = Nfolds, list = TRUE, returnTrain = FALSE)
fold_samples <- list()
for (k in 1:Nfolds) {
  fold <- flds[[k]]
  trains <- setdiff(1:N, fold)
  
  fit <- sampling(sm, data = list(N=N,d=d_obs,p=d_lat,Y=Y, seq=inseq, 
                                  Ntrain=length(trains), trains=trains, tests=fold), 
                  iter = 1, chains = 1)
  
  # Run GHMC
  max_it    <- 2000
  L <- 30
  
  fold_samples[[k]] <- run_GHMC(fit=fit, max_it=max_it, 
                      epsilons=list(sigma=.03, X=.008, u=.03, beta=.1, default.=0.01),
                      L=L, ignore_params=c("X", "u"), use_stiefel=TRUE,
                      N_generated_quantities = 2, thin=5, burnin=floor(max_it/2))
  
}

save(fold_samples, file=file.path("saved_samples", "rat_seq_crossval_samples.Rdata"))

# Compute predictions on test data
Nwrong = 0
for (k in 1:Nfolds) {
  fold <- flds[[k]]
  seq_true <- inseq[fold]
  print(seq_true)
  seq_probs <- laply(fold_samples[[k]], .fun=function(x) x$seq_probs)
  mean_probs <- colMeans(seq_probs)
  preds <- (mean_probs>.5)*1
  print(preds)
  print(seq_true == preds)
  Nwrong <- Nwrong + sum(seq_true !=preds)
}
print(Nwrong / N)


# compute expected log predictive densities on held out test data
Nsamples <- length(fold_samples[[1]])
elpds <- matrix(0., Nsamples, N)
for (k in 1:Nfolds) {
  fold <- flds[[k]]
  elpds[,fold] <- laply(fold_samples[[k]], .fun=function(x) x$lpd)
}
