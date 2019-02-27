set.seed(1)
dat <- readMat(file.path("data", "bit_vector.mat"))

true_bits <- dat$cleanX
corrupted <- dat$corrX
miss <- dat$newmiss
Y <- dat$X

N     <- dim(Y)[1]
d_obs <- dim(Y)[2]
d_lat <- 3

# compile Stan model for AD
sm  <- stan_model(file.path("stan", "epca.stan"))

fit <- sampling(sm, data = list(N=N,d=d_obs,p=d_lat,Y=Y, Ytrue=corrupted,miss=miss), 
              iter = 1, chains = 1)

# Run GHMC
epsilons=list(X=.03, u=.03, mu=0.03, lambda=.01)
L <- 80
max_it = 10000
joint_samples <- run_GHMC(fit=fit, max_it=max_it, 
                    epsilons=epsilons,
                    L=L, ignore_params=c("u"), use_stiefel=TRUE,
                    N_generated_quantities = 3, thin=1, burnin=0)

filen=file.path("saved_samples", sprintf("bit_vector_%d_samples.Rdata", max_it))
save(joint_samples, file=filen)
load(filen)

library(loo)
loglik <- laply(joint_samples, .fun = function(x) x$log_lik_miss)
preds <- laply(joint_samples, .fun = function(x) x$predicts)
alllik <- laply(joint_samples, .fun = function(x) x$alllik)
Xsamples <- laply(joint_samples, .fun = function(x) x$X)
Nsamples <- dim(Xsamples)[1]

predicts <- (exp(preds) > 0.5)
vs_true <- predicts
vs_corrupt <- predicts
for (s in 1:dim(preds)[1]) {
  vs_true[s,,] = 1*(predicts[s,,] != true_bits)
  vs_corrupt[s,,] = 1*(predicts[s,,] != corrupted)
}
vs_true <- rowMeans(vs_true)
vs_corrupt <- rowMeans(vs_corrupt)
par(mfrow=c(1,2))
plot(vs_true)
plot(vs_corrupt)

reconstructed <- 1 * (colMeans(exp(preds[100:500,,]))>0.5)
print(mean(reconstructed != true_bits))
reconstructed <- 1 * (colMeans(exp(preds[100:Nsamples,,]))>0.5)
print(mean(reconstructed != true_bits))

res <- list(mean_err=mean(reconstructed != true_bits),
            reconstruction=reconstructed)


X_0 <- chordal_mean(Xsamples[(Nsamples/2): Nsamples,,])
dists <- array(0, Nsamples)
for (j in 1:Nsamples)
  dists[j] <- dist_chordal(Xsamples[j,,], X_0)
#plot(dists)

#reconstruction_errors <- vapply(1:20, get_error, FUN.VALUE=1)
save(true_bits, corrupted, res, Y, dists, vs_true, vs_corrupt, 
	file=file.path("saved_samples", "bit_vec_Stiefel_results_3.Rdata"))


