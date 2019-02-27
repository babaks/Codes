library(rstan)
library(MASS)
library(rstiefel)
library(expm)

SuperChris1 <- read.csv(file.path("data", "counts_superchris_1.csv"), header=FALSE)
Y <- as.matrix(SuperChris1)
Y <- Y

inseq <- read.csv(file.path("data", "seq_superchris_1.csv"), header=FALSE)
inseq   <- as.matrix(inseq)
colnames(inseq) <- c()
inseq   <- as.vector(inseq)
N     <- dim(Y)[1]
d_obs <- dim(Y)[2]
d_lat <- 5
K     <- 2

# compile Stan model for AD
seq_sm  <- stan_model(file.path("stan", "joint_rat_logit.stan"))
fit_seq <- sampling(seq_sm, data = list(N=N,d=d_obs,p=d_lat,Y=Y, odr=inseq), 
                iter = 1, chains = 1)

attributes(fit_seq)$inits[[1]]$X <- rmf.matrix(attributes(fit_seq)$inits[[1]]$X)
attributes(fit_seq)$inits[[1]]$sigma <- 0.01 * attributes(fit_seq)$inits[[1]]$sigma 

# Run GHMC
max_it    <- 10000
L <- 30

joint_samples <- run_GHMC(fit=fit_seq, max_it=max_it, 
                    epsilons=list(sigma=.03, X=.008, u=.03, beta=.1, beta_0=.1, default.=0.01),
                    L=L, ignore_params=c(), use_stiefel=TRUE,
                    N_generated_quantities = 0, thin=10, burnin=max_it/2)


save(joint_samples, file=file.path("saved_samples", "thinned_joint_samples_seq.Rdata"))


