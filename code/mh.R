# Metropolis Hastings -------------------------------------------------------------------------
library(tibble)
library(ggplot2)

# helpers -----------------------------------------------------------------
pdf_mvnorm <- function(x, mu, Sigma, P = NULL, d) {
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1 / D
    sqrtP <- D.inv * tmp$vt
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  K <- (2*pi)**-(d/2)
  p <- det(Sigma)^-0.5  *exp(-0.5*crossprod(sqrtP %*% (x - mu)))
  return(p*K)
}

pdf_mvstudent <- function(x, nu, mu, Sigma, P = NULL, d) {
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1 / D
    sqrtP <- D.inv * tmp$vt
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  K <- gamma(0.5*(nu+d))/(gamma(0.5*nu)*(nu*pi)**(0.5*d))
  p <- det(Sigma)**-0.5*(1+nu**-1*crossprod(sqrtP %*% (x - mu)))**(-0.5*(nu+d))
  return(p*K)
  }

draw_mvn <- function(mu, Sigma, d, P=NULL){
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    sqrtSigma <- D*tmp$vt
  }
  
  x <- mu + t(sqrtSigma)%*%matrix(rnorm(d))
  return(x)
}

# Metropolis Hastings -------------------------------------------------------------------------
metropolis_hastings <- function(target,
                                proposal,
                                draw_proposal,
                                niter,
                                d) {
  theta <- matrix(NA, nrow = d, ncol = niter)
  alpha <- vector(mode = "numeric", length = niter)
  
  # 1.) draw initial value
  # 2.) repeat fot i in 1,...,niter
  theta_last <- matrix(c(1,1), nrow=2)
  
  for (i in 1:niter) {
    
    # a.) draw caindidate theta*
    theta_cand <- draw_proposal()
    
    # b.) calc alpha
    log_alpha <- log(target(theta_cand)) +
      log(proposal(theta_last)) -
      log(target(theta_last)) -
      log(proposal(theta_cand))
    
    alpha[i] <- exp(log_alpha)
    
    # c.)
    #browser()
    if (alpha[i] >= 1) {
      # accept draw
      theta[, i] <- theta_cand
      theta_last <- theta_cand
    } else {
      if (alpha[i] > runif(1)) {
        # accept with prop alpha
        theta[, i] <- theta_cand
        theta_last <- theta_cand
      }
      else {
        # reject
        theta[, i] <- theta_last
      }
    }
  }
  
  return(list(theta = theta, alpha = alpha))
}


# target 
target <-function(x){
  p <- pdf_mvstudent(x,
              nu = 2,
              mu = matrix(c(0, 0)),
              Sigma = matrix(c(4.59, 3.88, 3.88, 4.59), nrow=2),
              d = 2)
  return(p)
  }

# proposal
proposal <- function(x){
  p <- pdf_mvnorm(x,
                  mu = matrix(c(0, 0)),
                  Sigma = diag(2),
                  d = 2,
                  P=NULL)
  return(p)
  }

# draw from proposal
draw_proposal <- function(){
  x <- draw_mvn(mu = matrix(c(0, 0)),
                Sigma = diag(2),
                d = 2,
                P=NULL)
  return(x)
}

# sample
resu <- metropolis_hastings(target, proposal, draw_proposal, niter=10000, d=2)

# MCMC estimates
rowMeans(resu$theta)
var(t(resu$theta))

mean(pmin(1, resu$alpha))

plot(t(resu$theta))



# plot bivariate mvn
e1 <- matrix(c(6 , 6))
e2 <- matrix(c(0.5, -0.5))
crossprod(e1, e2)
v1 <- sqrt(crossprod(e1)[1])
v2 <- sqrt(crossprod(e2)[1])
e1 <- e1/v1
e2 <- e2/v2
Q <- cbind(e1, e2)
Lambda <- diag(c(v1, v2))
Sigma <- Q%*%Lambda%*%t(Q)

# L <- matrix(c(1, -0.9, 0, -1), nrow = 2, ncol = 2)
# Sigma <- tcrossprod(L)
# Sigma <- diag(2)
# Sigma <- matrix(c(4.0, -9.6, -9.6, 64.0), nrow=2)
mu <- matrix(c(0,0))

nsim <- 20000
x <- matrix(0, nrow=2, ncol=nsim)
for (i in 1:nsim){
  x[,i] <- draw_mvn(mu, Sigma, 2)
}
cov(t(x))


data <- tibble(x1 = x[1, ], x2=x[2, ])

ggplot(data, aes(x=x1, y=x2))+
  geom_bin2d()

ggplot(data, aes(x=x1, y=x2))+
  geom_density2d()

ggplot(data, aes(x=x1, y=x2))+
  geom_hex()
