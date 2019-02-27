run_GHMC <- function (fit, max_it=500, epsilons=list(), L=20, ignore_params=NULL, use_stiefel="Stiefel", 
                      N_generated_quantities=0, thin=1, burnin=0) {
  # Code to initialize and GHMC sample from stan model.
  # Returns list of samples
  #
  # NOTE: fit MUST include the parameter named "X", 
  # the grassmann matrix, as it's final parameter.
  # NOTE2: fit CANNOT have any transformed parameters.
  #
  # fit:      compiled stan model. 
  # max_it:   number of iterations to run
  # epsilons: list of (param_name: epsilon) step sizes for GHMC to take 
  #     for each parameter. Can have 'default.' entry to set generic epsilon.
  # L:        number of steps per GHMC iteration
  # ignore_params: which parameters to not save samples for.
  # use_stiefel:  if TRUE, use Stiefel manifold sampler instead of Grassmann sampler
  # N_generated_quantities = number of additional parameters after X to remove
  # from the GHMC sampling steps
  # thin = k: only keep every k'th sample
  
  
  # Grab all model parameter info from "fit"
  fit_attr <- attributes(fit)
  
  param_names <- fit_attr$model_pars
  param_names <- param_names[1:(length(param_names)-1)] # remove lp__
  param_inits <- fit_attr$inits[[1]]
  if ("sigma" %in% param_names) {
    print("adjusting sigma parameter")
    param_inits$sigma <- param_inits$sigma * 0.1
  }
  # initialize X to a random orthonormal matrix
  d <- diag(t(param_inits$X) %*% param_inits$X)
  if (any(abs(d - 1) > 0.000001)){
    print("initial X is not an orthogonal matrix -- resetting a random initial value")
    param_inits$X <- rmf.matrix(param_inits$X)
  }
  
  # create list of indices for unconstrained params
  lengths <- vapply(param_inits[1:(length(param_inits)-N_generated_quantities)], length, FUN.VALUE=1)
  ind_high <- cumsum(lengths)
  ind_low <- ind_high - lengths + 1
  ind = list()
  for (pname in param_names[1:(length(param_names) - N_generated_quantities)]) {
    ind[[pname]] <- ind_low[[pname]]:ind_high[[pname]]
  }
  param_inits <- constrain_pars(fit, unconstrain_pars(fit, param_inits))
  
  # Define epsilons
  if ('default.' %in% names(epsilons)){
    default. <- epsilons$default.
    epsilons$default. <- NULL
  } else {
    default. <- 0.01
  }
  eps <- as.list(rep(default., length(param_names)))
  names(eps) <- param_names
  
  for (pname in intersect(names(epsilons), param_names)) {
    eps[[pname]] <- epsilons[[pname]]
  }
  remaining <- setdiff(names(epsilons), names(eps))
  if (length(remaining) > 0) {
    warning('parameters ', remaining, 
            ' were not found in the stan model, but were declared in epsilons.\n',
            immediate. = TRUE)
  }
  
  # Init sample storage as list (sample number) of lists (param name) of param values
  keep_pars <- setdiff(param_names, ignore_params)  # don't save unwanted params
  samples <- rep_len(list(param_inits[keep_pars]), floor((max_it-burnin) / thin))
  
  accept    <- rep(0, max_it) #keep track of acceptance
  auto_dist <- rep(0, max_it) #geodesic distance between samples
  
  #begin HMC
  blocksize <- 20
  param_prev <- param_inits
  if (max_it > 1) {
    i=1
    while (i < max_it){
      i = i + 1
      res <- try(doGHMC(theta = param_prev, eps=eps, L=L, fit=fit, ind=ind, 
                    use_stiefel=use_stiefel))
      
      if (!inherits(res, "try-error")) {
        accept[i] <- res$accept
        
        if (i %% thin == 0 && i > burnin)
          samples[[floor((i-burnin) / thin)]] <- res$theta[keep_pars]
        param_prev <- res$theta
        
        if (i %% blocksize ==0) print(c(i, mean(accept[((i-blocksize)+1):i])))
      } else { # Error in sampling so repeat previous sample
        if (i %% thin == 0 && i > burnin)
          samples[[floor((i-burnin) / thin)]] <- param_prev
      }
      
      if (i==2 && accept[i] == 0) { # Reset starting values if can't take first step.
        i = 1
        print("retrying first iteration")
        nunc <- get_num_upars(fit)
        param_prev <- constrain_pars(fit, .1 * runif(nunc) + 0.00001)
        param_prev$X <- rmf.matrix(attributes(fit)$inits[[1]]$X)
      }
    }
    
    print(sprintf("accept. rate: %0.3f", mean(accept)))
  }
  
  return(samples)
  
}

dist_grass <- function(X,Y) {
  svd.obj <- svd(t(X)%*%Y)
  d       <- svd.obj$d
  d[which(d>=1)] <- 1
  thetas    <- acos(d)
  return(sqrt(sum(thetas^2)))
}

chordal_mean <- function(Xs) {
  N <- dim(Xs)[1]
  Xhat <- Xs[1,,] %*% t(Xs[1,,])
  for (j in 2:N)
    Xhat <- Xhat + Xs[j,,] %*% t(Xs[j,,])
  U <- (svd(Xhat)$u)[, 1:(dim(Xs)[3])]
  return(U)
}

dist_chordal <- function(X, Y) {
  svd.obj <- svd(t(X)%*%Y)
  d       <- svd.obj$d
  d[which(d>=1)] <- 1
  thetas    <- acos(d)
  return(sqrt(sum(sin(thetas)^2)))
}

dist_to_chordal_mean <- function(samples, X0=NULL) {
  
  Nsamp <- length(samples)
  Xs <- laply(samples, .fun=function(x) x$X)
  if (is.null(X0))
    chmean <- chordal_mean(Xs[(Nsamp / 2):Nsamp,,])
  else
    chmean <- X0
  return(aaply(Xs, .fun=function(x) dist_chordal(x, chmean), .margins=1))
}

cnpp1 <- function(n,p) {
  res <- 1/gamma(p*(n-p)/2 + 1)
  for (j in 1:p) {
    res <- res * gamma((n-j+1)/2) / gamma((p-j+1)/2)
  }
  return(res)
}

upper_bnd_mu <- function(delta, n, p) {
  if (abs(delta) > 1)
    return(1)
  
  return(cnpp1(n,p) * delta^(p*(n-p)) / (1-delta^2)^(p/2))
}

