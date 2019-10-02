library(matrixcalc)
library(mvtnorm)

# Metropolis Hastings -------------------------------------------------------------------------


NIG_prior <- function(theta, mu_0, sigma_0, nu_0, beta_0) {

  # Normal prior for mean
  N <- pnorm(x, nu_0, sigma_0)

  # Inverse Gamma prior for variance
  IG <- 1 / dgamma(x, nu_0, beta_0)

  return(N * IG)
}

likelihood <- function(X, mu, sigma) {

  # likelihood
  l <- dmvnorm(X, mu, sigma)

  return(cumprod(l))
}


# posterior <- function() {
#
#
#
# }

trunc_norm <- function(x) {
  if (!is.null(dim(x))) {
    d <- max(dim(x))
  }
  else {
    d <- length(x)
  }

  return(
    dmvnorm(abs(x), mean = matrix(0, nrow = d, ncol = 1), sigma = diag(d))
  )
}

n <- 10000
A <- cbind(seq(-2, 2, 1e-3), seq(-2, 2, 1e-3))
resu <- vector(mode = "numeric", length = nrow(A))
for (i in 1:nrow(A)) {
  resu[i] <- trunc_norm(A[i, ])
}


target <- function(x) {
  return(dbeta(x, shape1 = 2, shape2 = 5))
}

proposal <- function(x=0, mean=0, type="pdf") {
  # rmvnorm(1, mean = matrix(0, nrow = d, ncol = 1))
  
  if (type=="pdf") {
    return(1)
  } else if (type=="rv"){
    return(runif(1))
  }
  
  
}




# Metropolis-Hastings
metropolis_hastings <- function(target, proposal, niter = 100, d) {
  theta <- matrix(NA, nrow = niter, ncol = d)
  alpha <- vector(mode="numeric", length = niter)
  # 1.) draw initial value
  # 2.) repeat fot i in 1,...,niter
  theta_last <-  proposal(x=0, mean = theta_last, type="rv")
    for (i in 1:niter) {

      # a.) draw caindidate theta*
      theta_cand <- proposal(x=0, mean = theta_last, type="rv")
      # dmvnorm(theta_last, mean=theta_cand)
      # dmvnorm(theta_cand, mean=theta_last)
      # b.) calc alpha
      alpha[i] <- (target(theta_cand) * proposal(theta_last, mean = theta_cand)) / (target(theta_last) * proposal(theta_cand, mean = theta_last))
      # c.)

      if (alpha[i] >= 1) {
        # accept draw
        theta[i] <- theta_cand
        theta_last <- theta_cand
      } else {
        if (alpha[i] > runif(1)) {
          theta[i] <- theta_cand
          theta_last <- theta_cand
        }
        else {
          # accept with prop alpha
          theta[i] <- theta_last
          # theta_last=theta_cand
        }
      }
    }
  
  return(list(theta = theta, alpha = alpha))
}


mh <- metropolis_hastings(target, proposal, niter = 1e4, d=1)

hist(mh$theta)
mean(mh$theta)
var(mh$theta)


mean_beta <- cumsum(mh$theta)/c(1:length(mh$theta))

plot.ts(mean_beta)



