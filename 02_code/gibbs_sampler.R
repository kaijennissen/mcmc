library(matrixcalc)
library(mvtnorm)

# Metropolis Hastings -------------------------------------------------------------------------
target <- function(x) {
  return(dbeta(x, shape1 = 2, shape2 = 5))
}



mh <- metropolis_hastings(target, proposal, niter = 1e4, d=1)

hist(mh$theta)
mean(mh$theta)
var(mh$theta)


mean_beta <- cumsum(mh$theta)/c(1:length(mh$theta))

plot.ts(mean_beta)



# NIG Example
beta <- c(1, 0.5, -3)
sigma <- 10
X <- matrix(runif(3000, -100, 100), nrow=1000)
y <- X%*%beta + rnorm(1000, sd = sigma)


gamma_prior <- list(a0 = 1e-3, b0 = 1e-4)
beta_prior <- list(beta0 = matrix(0, nrow=3, ncol=1), sigma0 = 300)
 
# Gibbs-Sampler
gibbs_NIG <- function(y, X ,nsim, beta_prior, gamme_prior){

  beta0 <- beta_prior$beta0
  sigma2_0 <- beta_prior$sigma0
  a0 <- gamma_prior$a0
  b0 <- gamma_prior$b0
  N <- length(y)
  k <- dim(beta0)[1]
  I3 <- diag(k)
  gamma <- 15
  new_a <- a0+0.5*crossprod(y-X%*%beta)
  sigma2 <- 3000
  
  beta_store <- matrix(0, k, nsim)
  sigma2_store <- vector("numeric", nsim)
  
  for (i in 1:nsim){
  
  Sigma0_inv <- I3/sigma2_0
  Sigma_inv <- diag(N)/sigma2
  Sigma_new <- (crossprod(X, Sigma_inv)%*%X + Sigma0_inv)
  beta_hat <- solve(Sigma_new, crossprod(X, Sigma_inv)%*%y+Sigma0_inv%*%beta0)
  L <- chol(Sigma_new)
  
  beta<- beta_hat + solve(L, rnorm(k))
  
  new_b <- b0+0.5*crossprod(y-X%*%beta)
  sigma2 <- rgamma(1, shape = new_a, rate = new_b)[]
  
  beta_store[,i] <- beta
  sigma2_store[i] <- sigma2
  i=i+1
  }
  
  return(list(beta_store, sigma2_store))
}

resu <- gibbs_NIG(y, X, nsim=1000, beta_prior, gamme_prior)


beta <- resu[[1]]
sigma2 <- resu[[2]]

rowMeans(beta)
mean(sigma2)


hist(beta[1,])
hist(beta[2,])
hist(beta[3,])
hist(sigma2)

plot(beta[1,], beta[2,])
plot(beta[1,], beta[3,])
plot(beta[2,], beta[3,])
plot(beta[1,], sigma2)
plot(beta[2,], sigma2)
plot(beta[3,],sigma2)



# MH Sampler --------------------------------------------------------------
log_target <- function(y, X, theta, beta0, Sigma0, a0, b0){
  beta <- theta[1:3, ,FALSE]
  sigma2 <- theta[4, ,TRUE]
  
  N <- length(y)
  p <- log_mvnorm(y, X%*%beta, sigma2*diag(N))+
    log_mvnorm(beta, beta0, Sigma0)+
    log_gamma(sigma2, a0, b0)
  return(p)
}

log_mvnorm <- function(x, mu, Sigma, P=NULL){
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1/D
    sqrtP <- D.inv * tmp$vt 
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  p <- -0.5*log(det(Sigma)) + crossprod(sqrtP %*% (x-mu))
  return(p)
}

log_gamma <- function(x, a, b){
  p <- (a-1)*x-x*b
  return(p)
}

log_proposal <- function(x=0, mu=0, Sigma=0) {
  p <- log_mvnorm(x, mean, Sigma)
  return(p)
}

draw_proposal <- function(x=0, mu=0, Sigma=0, P=NULL) {
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1/D
    sqrtP <- D.inv * tmp$vt 
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  k <- dim(mu)[1]
  x <- mu+sqrtP%*%rnorm(k)
  return(x)
}

# Metropolis-Hastings
metropolis_hastings <- function(y, X
                                log_target,
                                log_proposal,
                                draw_proposal,
                                niter = 100,
                                d,
                                beta0,
                                Sigma0, 
                                a0,
                                b0) {
  
  theta <- matrix(NA, nrow = niter, ncol = d)
  alpha <- vector(mode="numeric", length = niter)
  
  # 1.) draw initial value
  # 2.) repeat fot i in 1,...,niter
  theta_last <-  draw_proposal(mean = matrix(0, ncol=d), 1e5*diag(d))
  
  for (i in 1:niter) {
    
    # a.) draw caindidate theta*
    theta_cand <- draw_proposal(mean = theta_last, 1e5*diag(d))
    
    # b.) calc alpha
    function(y, X, theta_cand, beta0, Sigma0, a0, b0)
      log_alpha <-  log_target(y, X, theta_cand, beta0, Sigma0) +
        log_proposal(theta_last, mean = theta_cand) -
        log_target(y, X, theta_last, beta0, Sigma0, a0, b0) -
        log_proposal(theta_cand, mean = theta_last)
      alpha[i] <-exp(log_alpha)
      
      # c.)
      if (alpha[i] >= 1) {
        # accept draw
        theta[,i] <- theta_cand
        theta_last <- theta_cand
      } else {
        if (alpha[i] > runif(1)) {
          # accept with prop alpha
          theta[,i] <- theta_cand
          theta_last <- theta_cand
        }
        else {
          # reject
          theta[,i] <- theta_last
        }

        
      }
  }
  
  return(list(theta = theta, alpha = alpha))
}


function(y, X
         log_target,
         log_proposal,
         draw_proposal,
         niter = 100,
         d,
         beta0,
         Sigma0, 
         a0,
         b0)



