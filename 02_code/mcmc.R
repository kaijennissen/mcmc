# Metropolis Hastings ---------------------------------------------------------
# beta ditribution
target <- function(x) {
  return(dbeta(x, shape1 = 2, shape2 = 5))
}

# Uniform distribution
proposal <- function(x = 0, mean = 0, type = "pdf") {
  
  if (type == "pdf") {
    return(1)
  } else if (type == "rv"){
    return(runif(1))
  }
  
}


# Metropolis-Hastings Sampler
MH_sampler <- function(target, proposal, niter = 100, d) {
  
  # storage
  theta <- matrix(NA, nrow = niter, ncol = d)
  alpha <- vector(mode="numeric", length = niter)
 
   # 1.) draw initial value
  # 2.) repeat fot i in 1,...,niter
  theta_last <-  proposal(x = 0, mean = theta_last, type = "rv")
    
  for (i in 1:niter) {

      # a.) draw caindidate theta*
      theta_cand <- proposal(x=0, mean = theta_last, type = "rv")
     
       # b.) calc alpha
      log_alpha = (log(target(theta_cand)) + 
        log(proposal(theta_last, mean = theta_cand)) - 
          log(target(theta_last)) -
          log(proposal(theta_cand, mean = theta_last))
      )     
       alpha[i] <- exp(log_alpha)
      
      # c.)
      if (alpha[i] >= 1) {
        # accept draw wit prop = 1
        theta[i] <- theta_cand
        theta_last <- theta_cand
      } else {
        # accept draw with prop = alpha
        if (alpha[i] > runif(1)) {
          theta[i] <- theta_cand
          theta_last <- theta_cand
        }
        else {
          theta[i] <- theta_last
        }
      }
    }
  
  return(list(theta = theta, alpha = alpha))
}


# run MH-Sampler
mh <- MH_sampler(target, proposal, niter = 1e5, d=1)



hist(mh$theta)
mean(mh$theta) # true mean: 0.2857143
var(mh$theta) # # true variance: 0.0255102


mean_beta <- cumsum(mh$theta)/c(1:length(mh$theta))

plot.ts(mean_beta)
plot.ts(cumsum(pmin(mh$alpha, 1))/c(1:length(mh$alpha)))

       