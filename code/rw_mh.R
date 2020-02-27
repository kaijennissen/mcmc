# Metropolis Hastings -------------------------------------------------------------------------
library(tibble)
library(ggplot2)

# helpers -----------------------------------------------------------------
log_pdf_mvnorm <- function(x, mu, Sigma, P = NULL, d) {
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1 / D
    sqrtP <- D.inv * tmp$vt
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  K <- -(d/2)*log((2*pi))
  p <- -0.5*log(det(Sigma))-0.5*crossprod(sqrtP %*% (x - mu))
  return(p+K)
}

log_pdf_mvstudent <- function(x, nu, mu, Sigma, P = NULL, d) {
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1 / D
    sqrtP <- D.inv * tmp$vt
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  K <- log(gamma(0.5*(nu+d)))-log(gamma(0.5*nu))+(0.5*d)*log((nu*pi))
  p <- -0.5*log(det(Sigma))+(-0.5*(nu+d)*log(1+nu**-1*crossprod(sqrtP %*% (x - mu))))
  return(p+K)
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
random_walk_metropolis_hastings <- function(target,
                                proposal,
                                draw_mvn,
                                niter,
                                d) {
  theta <- matrix(NA, nrow = d, ncol = niter)
  alpha <- vector(mode = "numeric", length = niter)
  
  # 1.) draw initial value
  # 2.) repeat fot i in 1,...,niter
  theta_curr <- matrix(c(0.5, 0.5), nrow=2)
  
  for (i in 1:niter) {
    
    # a.) draw caindidate theta*
    theta_prop <- draw_mvn(mu = theta_curr, Sigma = 3*diag(2), d=2)
    
    # b.) calc alpha
    log_alpha <- log_target(x = theta_prop) -
      log_target(x = theta_curr) 
    
    alpha[i] <- exp(log_alpha)
    
    # c.)
    #browser()
    if (alpha[i] >= 1) {
      # accept draw
      theta[, i] <- theta_prop
      theta_curr <- theta_prop
    } else {
      if (alpha[i] > runif(1)) {
        # accept with prop alpha
        theta[, i] <- theta_prop
        theta_curr <- theta_prop
      }
      else {
        # reject
        theta[, i] <- theta_curr
      }
    }
  }
  
  return(list(theta = theta, alpha = alpha))
}


# target 
log_target <-function(x){
  p <- log_pdf_mvstudent(x,
              nu = 12,
              mu = matrix(c(0, 0)),
              Sigma = matrix(c(4.59, 3.88, 3.88, 4.59), nrow=2),
              d = 2)
  return(p)
  }

# proposal
log_proposal <- function(x){
  p <- log_pdf_mvnorm(x,
                  mu = matrix(c(0, 0)),
                  Sigma = 3*diag(2),
                  d = 2,
                  P=NULL)
  return(p)
  }


# sample
resu <- random_walk_metropolis_hastings(target, proposal, draw_mvn, niter=3e5, d=2)
nu <- 12
# MCMC estimates
# mean  = mu
# var = nu/(nu-2)*Sigma
rowMeans(resu$theta)
var(t(resu$theta))
nu/(nu-2)*Sigma

mean(pmin(1, resu$alpha))

draw_mh <- tibble(dist = "rw_mh", x1 = resu$theta[1,], x2 = resu$theta[2,])





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
# Sigma <-matrix(c(4.59, 3.88, 3.88, 4.59), nrow=2)# matrix(c(4.0, -9.6, -9.6, 64.0), nrow=2)
mu <- matrix(c(0,0))

nsim <- 1e6
x <- matrix(0, nrow=2, ncol=nsim)
for (i in 1:nsim){
  x[,i] <- draw_mvn(mu, Sigma, 2)
}
cov(t(x))


data_norm <- tibble(dist = "norm", x1 = x[1, ], x2=x[2, ])

# ggplot(data_norm, aes(x=x1, y=x2))+
#   geom_bin2d()
# 
# ggplot(data_norm, aes(x=x1, y=x2))+
#   geom_density2d()
# 
# ggplot(data_norm, aes(x=x1, y=x2))+
#   geom_hex()

# multivariate t
draw_mvt <- function(n, nu, mu, Sigma){
  d = dim(mu)[1]
  x <- matrix(0, nrow=2, ncol=n )
  for (i in 1:n){
    z <- sum(rnorm(nu)^2)
    x[,i] <- draw_mvn(mu, Sigma, d=2) * sqrt(nu / z) + mu
    
  }
  return(t(x))
}

Sigma <- matrix(c(4.59, 3.88, 3.88, 4.59), nrow=2)
x <- draw_mvt(1e5, 12, matrix(0, 2), Sigma)

data_t <- tibble(dist = "t", x1 = x[,1], x2=x[,2])   

# ggplot(data_t, aes(x=x1, y=x2))+
#   geom_bin2d()   
# 
# ggplot(data_t, aes(x=x1, y=x2))+
#   geom_density2d()
# 
# ggplot(data_t, aes(x=x1, y=x2))+
#   geom_hex()


data <- dplyr::bind_rows(draw_mh, data_t)

ggplot(data[abs(data$x1)<12 & abs(data$x2)<12, ], aes(x=x1, y=x2))+
  geom_hex()+
  facet_wrap(~dist)



   