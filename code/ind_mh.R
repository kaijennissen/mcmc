# We'll assume estimation of a Poisson mean as a function of x
x <- runif(100)
y <- rpois(100,5*x)  # beta = 5 where mean(y[i]) = beta*x[i]

# Prior distribution on log(beta): t(5) with mean 2 
# (Very spread out on original scale; median = 7.4, roughly)
log_prior <- function(log_beta) dt(log_beta-2, 5, log=TRUE)

# Log likelihood
log_lik <- function(log_beta, y, x) sum(dpois(y, exp(log_beta)*x, log=TRUE))

# Random Walk Metropolis-Hastings 
# Proposal is centered at the current value of the parameter

rw_proposal <- function(current) rnorm(1, current, 0.25)
rw_p_proposal_given_current <- function(proposal, current) dnorm(proposal, current, 0.25, log=TRUE)
rw_p_current_given_proposal <- function(current, proposal) dnorm(current, proposal, 0.25, log=TRUE)

rw_alpha <- function(proposal, current) {
	# Due to the structure of the rw proposal distribution, the rw_p_proposal_given_current and
	# rw_p_current_given_proposal terms cancel out, so we don't need to include them - although
	# logically they are still there:  p(prop|curr) = p(curr|prop) for all curr, prop
	exp(log_lik(proposal, y, x) + log_prior(proposal) - log_lik(current, y, x) - log_prior(current))
}

# Independent Metropolis-Hastings
# Note: the proposal is independent of the current value (hence the name), but I maintain the
# parameterization of the functions anyway.  The proposal is not ignorable any more
# when calculation the acceptance probability, as p(curr|prop) != p(prop|curr) in general.

ind_proposal <- function(current) rnorm(1, 2, 1) 
ind_p_proposal_given_current <- function(proposal, current) dnorm(proposal, 2, 1, log=TRUE)
ind_p_current_given_proposal <- function(current, proposal) dnorm(current, 2, 1, log=TRUE)

ind_alpha <- function(proposal, current) {
	exp(log_lik(proposal, y, x)  + log_prior(proposal) + ind_p_current_given_proposal(current, proposal) 
		- log_lik(current, y, x) - log_prior(current) - ind_p_proposal_given_current(proposal, current))
}

# Vanilla Metropolis-Hastings - the independence sampler would do here, but I'll add something
# else for the proposal distribution; a Normal(current, 0.1+abs(current)/5) - symmetric but with a different
# scale depending upon location, so can't ignore the proposal distribution when calculating alpha as
# p(prop|curr) != p(curr|prop) in general

van_proposal <- function(current) rnorm(1, current, 0.1+abs(current)/5)
van_p_proposal_given_current <- function(proposal, current) dnorm(proposal, current, 0.1+abs(current)/5, log=TRUE)
van_p_current_given_proposal <- function(current, proposal) dnorm(current, proposal, 0.1+abs(proposal)/5, log=TRUE)

van_alpha <- function(proposal, current) {
	exp(log_lik(proposal, y, x)  + log_prior(proposal) + ind_p_current_given_proposal(current, proposal) 
		- log_lik(current, y, x) - log_prior(current) - ind_p_proposal_given_current(proposal, current))
}


# Generate the chain
values <- rep(0, 10000) 
u <- runif(length(values))
naccept <- 0
current <- 1  # Initial value
propfunc <- van_proposal  # Substitute ind_proposal or rw_proposal here
alphafunc <- van_alpha    # Substitute ind_alpha or rw_alpha here
for (i in 1:length(values)) {
	proposal <- propfunc(current)
	alpha <- alphafunc(proposal, current)
	if (u[i] < alpha) {
		values[i] <- exp(proposal)
		current <- proposal
		naccept <- naccept + 1
	} else {
		values[i] <- exp(current)
	}
}
naccept / length(values)
summary(values)
