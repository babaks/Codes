# STIEFEL projection and flow

stiefel_proj <- function(X,V){ #proj V onto tangent space of Stiefel at X
  proj <- V - .5 *X%*%(t(X)%*%V + t(V)%*%X)
  return(proj)
}

stiefel_flow <- function(Y_0,Y_dot_0,t){ #geodesic flow for t time
  p       <- dim(Y_0)[2]
  X_0_V_0 <- cbind(Y_0,Y_dot_0)
  A       <- t(Y_0)%*%Y_dot_0
  M_top   <- cbind(A,-t(Y_dot_0)%*%Y_dot_0)
  M_bot   <- cbind(diag(p),A)
  M       <- t*rbind(M_top,M_bot)
  M_1     <- expm::expm(M, method="Pade", order=6)
  M       <- expm::expm(-t*A, method="Pade", order=6)
  zeros   <- matrix(0, dim(M)[1], dim(M)[2])
  M_2     <- rbind(cbind(M, zeros), cbind(zeros, M))
  M_3     <- as.matrix(X_0_V_0%*%M_1%*%M_2)
  Y_t     <- M_3[,1:p]
  
  # Re-orthogonalize Y_t
  sv <- svd(Y_t)
  Y_t <- sv$u %*% t(sv$v)
    
  Y_dot_t <- M_3[,(p+1):(2*p)]
  return(list(Y_t,Y_dot_t))
}

# GRASSMANN flow and projection

grass_flow <- function(Y_0,Y_dot_0,t){
  svd_Y_dot_0 <- svd(Y_dot_0)
  U     <- svd_Y_dot_0$u
  V     <- svd_Y_dot_0$v
  
  Sigma <- diag(svd_Y_dot_0$d)
  Y_t     <- Y_0 %*% V %*% diag(cos(diag(Sigma)*t)) %*% t(V) +
    U %*% diag(sin(diag(Sigma)*t)) %*% t(V)
  Y_dot_t <- Y_0 %*% V %*% Sigma %*% diag(-sin(diag(Sigma)*t)) %*% t(V) +
    U %*% Sigma %*% diag(cos(diag(Sigma)*t)) %*% t(V)
 
  sv <- svd(Y_t)
  Y_t <- sv$u %*% t(sv$v)
  Y_dot_t <- grass_proj(Y_t, Y_dot_t)
  
  distance <- sqrt(sum(diag(Sigma)^2)) * t
  
  return(list(Y_t,Y_dot_t,distance))
}

grass_proj <- function(X,V){
  out <- (diag(dim(V)[1]) - X%*%t(X))%*%V
  return(out)
}

# Euclidean flow and project

euclid_flow <- function(Y_0, Y_dot_0, t) {
  return(Y_0 + t * Y_dot_0)
}

euclid_proj <- function(X, V) { 
  return(V)  
}


doGHMC <- function(theta, eps, L, fit, ind, use_stiefel="Stiefel") {
  
  if (use_stiefel==TRUE) {
    use_stiefel <- "Stiefel"
  } else if (use_stiefel ==FALSE) {
    use_stiefel <- "Grassmann"
  }
  if (!(use_stiefel %in% c("Euclid", "Stiefel", "Grassmann"))) {
    stop("select a valid sampler")
  }
  if (use_stiefel == 'Euclid') {
    flow    <- euclid_flow
    project <- euclid_proj
  } else if (use_stiefel == 'Stiefel') {
    flow    <- stiefel_flow
    project <- stiefel_proj
  } else {
    flow    <- grass_flow
    project <- grass_proj
  }
  
  tot_length <- get_num_upars(fit)
  length_but_X <- tot_length - length(theta$X)
  
  names_but_X <- setdiff(names(ind), 'X')
  
  # (1) draw V and p
  V <- matrix(rnorm(theta$X), dim(theta$X))
  p <- rnorm(length_but_X)
  
  # (2) project V onto tangent space at X
  V <- project(X=as.matrix(theta$X),V=V)
  V_vec <- as.vector(V)
  
  # (3) get log_prob
  un <- unconstrain_pars(fit, pars=theta) 
  h <- log_prob(fit,upars=un, adjust_transform=TRUE)[1] -
    .5*t(V_vec)%*%V_vec -.5*t(p)%*%p
  
  # (4)
  L <- round(L * (.1 + runif(1)))
  
  # (5)
  # (6) Half-step momentum for X
  grad  <- grad_log_prob(fit,upars=un,adjust_transform=TRUE)[1:tot_length]
  V_vec <- V_vec + eps$X / 2 * grad[ind$X]
  X_star <- matrix(un[ind$X], dim(theta$X))
  V     <- project(X=X_star, V=matrix(V_vec, dim(theta$X)))
  #browser()
  # Half-step momentum of remaining vars
  for (pname in names_but_X) {
    p[ind[[pname]]] <- p[ind[[pname]]] + eps[[pname]] / 2 * grad[ind[[pname]]]
  }
  
  for(tau in 1:L) {

    # (8) full step of geodesic flow for X
    X_t_V_t <- flow(Y_0=X_star,
                    Y_dot_0=V, eps$X)
    un[ind$X]  <- as.vector(X_t_V_t[[1]])
    V           <- X_t_V_t[[2]]
    V_vec       <- as.vector(V)
    
    # full step for position of remaining parameters
    for (pname in names_but_X) {
      un[ind[[pname]]] <- un[ind[[pname]]] + eps[[pname]] * p[ind[[pname]]]
    }
    
    # (9)
    if (tau < L) {
    # full step for momentum of remaining vars
      grad  <- grad_log_prob(fit,upars=un,adjust_transform=TRUE)[1:tot_length]
      V_vec <- V_vec + eps$X * grad[ind$X]
      X_star <- matrix(un[ind$X], dim(theta$X))
      V     <- project(X=X_star, V=matrix(V_vec, dim(theta$X)))
      V_vec <- as.vector(V)
      for (pname in names_but_X) {
        p[ind[[pname]]] <- p[ind[[pname]]] + eps[[pname]] * grad[ind[[pname]]]
      }
    }
  } # (11) end for
  
  # half step for momentum of remaining vars
  grad  <- grad_log_prob(fit,upars=un,adjust_transform=TRUE)[1:tot_length]
  V_vec <- V_vec + eps$X / 2 * grad[ind$X]
  X_star <- matrix(un[ind$X], dim(theta$X))
  V     <- project(X=X_star, V=matrix(V_vec, dim(theta$X)))
  V_vec <- as.vector(V)
  for (pname in names_but_X) {
    p[ind[[pname]]] <- p[ind[[pname]]] + eps[[pname]] / 2 * grad[ind[[pname]]]
  }
  
  theta_star <- constrain_pars(fit, upars=un)
  
  h_star <- log_prob(fit,upars=un, adjust_transform=TRUE)[1] -
    .5*t(V_vec)%*%V_vec -.5*t(p)%*%p
  
  # (13-16)
  z <- runif(1)
  if(z < exp(h_star-h)){
    return(list(accept=1, theta=theta_star))
  }else{
    return(list(accept=0, theta=theta))
  }
}
