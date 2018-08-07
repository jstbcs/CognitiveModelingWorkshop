
# Example of the Metropolis algorithm 
# Estimating the mean of a normal distribution with known variance (15^2)
# and a normal(mu_prior, sigma_prior) prior distribution on the mean

#### FUNCTIONS ----

metropolis_iq <- function(y, prior_mean, prior_sd, proposal_sd = 5, n_samples = 1000, start = 100){
  # metropolis algorithm for estimating the mean of a normal distribution
  # with a known variance (15 like iq tests) and a normal prior distribution on the mean
  
  samples = rep(NA, n_samples) # create a vector to store our samples
  current_mu = start # set the current value of mu to the starting value (set by user)
  
  for (i in 1:n_samples){ # do the samples
    
    samples[i] = current_mu # store the current value of mu
    
    proposed_mu = rnorm(n = 1, mean = current_mu, sd = proposal_sd) # sample a proposed value of mu
    # the proposal distribution is itself normal and is centered on the current value with a modifyable
    # standard deviation (which determines how much the chain can move around)
    
    # calculate the posterior probabilities for the current and proposed value of mu
    # we do things on the log scale to minimize problems due to rounding
    f_current = dnorm(x = current_mu, mean = prior_mean, sd = prior_sd, log = T) + sum(dnorm(x = y, mean = current_mu, sd = 15, log = T))
    f_proposed = dnorm(x = proposed_mu, mean = prior_mean, sd = prior_sd, log = T) + sum(dnorm(x = y, mean = proposed_mu, sd = 15, log = T))
    
    # take the ratio and thus the probability of accepting the proposed move
    # as we are working on the log scale we subtract rather then divide
    p_accept = exp(min(f_proposed - f_current, 0)) # same as min(f_proposed/f_current, 1) on natural scale
    
    if (runif(n = 1, min = 0, max = 1) < p_accept){
      # if the random uniform value is smaller than p_accept set the proposed mu 
      # to the current mu for the next step
      current_mu = proposed_mu
    }
  }
  # return the vector of sampled values of mu
  return(samples)
}

posterior_quantities <- function(y, prior_mean, prior_sd, known_sd = 15){
  # for a normal likelihood with a known variance and normal prior on the mean
  # the posterior distribution is also normal
  # this function returns the mean and sd of the posterior distribution
  # for comparison with the results of the metropolis algorithm
  
  # details on the derivation can be found here: https://mcs.utm.utoronto.ca/~nosedal/sta313/sta313-normal-mean.pdf
  ybar = mean(y)
  n = length(y)
  
  prior_prec = 1/prior_sd^2
  lik_prec = n/known_sd^2
  
  post_prec = prior_prec + lik_prec
  post_sd = sqrt(1/post_prec)
  
  post_mu = prior_mean*(prior_prec/post_prec) + ybar*(lik_prec/post_prec)
  
  return(c(post_mu, post_sd))
}


#### START HERE ----

N = 50 # number of observations in the group
mu_real = 110 # the 'true' mean

set.seed(123)
Y = rnorm(N, mu_real, 15) # simulate data
prior_m = 100
prior_sd = 10

# run the algorithm
samples = metropolis_iq(y = Y, # data
                        # set the first 2 arguments above
                        prior_mean = prior_m, # prior mean on the unknown mean (mu) we are trying to estimate
                        prior_sd = prior_sd, # prior standard deviation on mu
                        proposal_sd = 5, # standard deviation of the proposal distribution
                        n_samples = 1000, # number of samples from the posterior that we want
                        start = 100) # starting value for mu

# plot the steps in the chain
plot(samples, xlab = "Step", ylab = bquote(mu), type="l", col="blue")

# plot the histogram
hist(samples, xlab = bquote(mu), main = paste0("N samples = ", length(samples)), col = "grey", breaks = 30, probability = T)

# calculate the exact form of the posterior to compare to the samples
post = posterior_quantities(y = Y, prior_mean = prior_m, prior_sd = prior_sd)
curve(dnorm(x, mean = post[1], sd = post[2]), from = min(samples), to = max(samples), col = "red", lwd = 2, add = T)

# run more samples to better approximate the posterior

more_samples = metropolis_iq(y = Y, 
                        prior_mean = prior_m, 
                        prior_sd = prior_sd,
                        proposal_sd = 5, 
                        n_samples = 10000, # instead of 1000
                        start = 100)

plot(more_samples, xlab = "Step", ylab = bquote(mu), type="l", col="blue")

# plot the histogram
hist(more_samples, xlab = bquote(mu), main = paste0("N samples = ", length(more_samples)), col = "grey", breaks = 30, probability = T)
curve(dnorm(x, mean = post[1], sd = post[2]), from = min(samples), to = max(samples), col = "red", lwd = 2, add = T)

# an example showing the need for burn in or warm up

burn_samples = metropolis_iq(y = Y, 
                             prior_mean = prior_m, 
                             prior_sd = prior_sd,
                             proposal_sd = 5, 
                             n_samples = 2000,
                             start = 10)

# the starting value can influence the first steps of the chain
plot(burn_samples, xlab = "Step", ylab = bquote(mu), type="l", col="blue", main="before burn in")

# so was can disgard some samples as a burn in period
plot(x = 1001:2000, burn_samples[1001:2000], xlab = "Step", ylab = bquote(mu), type="l", col="blue", main = "after burn in")

# an example of a chain with bad autocorrelation

autocor_samples = metropolis_iq(y = Y, 
                             prior_mean = prior_m, 
                             prior_sd = prior_sd,
                             proposal_sd = .25, # this means that each step in the chain will be smaller. Leading to more autocorrelation
                             n_samples = 10000,
                             start = 100)

plot(autocor_samples, xlab = "Step", ylab = bquote(mu), type="l", col="blue", main="before thinning")

K = 10 # what level of thinning should we use? We'll keep every K^th sample
thin_samples = autocor_samples[seq(0, 10000, by = K)]

plot(thin_samples, xlab = "Step", ylab = bquote(mu), type="l", col="blue", main="after thinning")

coda::autocorr.plot(coda::as.mcmc(autocor_samples, main="before thinning"))

coda::autocorr.plot(coda::as.mcmc(thin_samples, main="after thinning"))

