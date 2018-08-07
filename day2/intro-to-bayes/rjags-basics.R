
library(rjags)
library(coda)

set.seed(123)
N = 100
mu_real = 110

Y = rnorm(N, mu_real, 15) # data

data_list = list('y' = Y,
            'N' = N)

model = "
  model {
    for (i in 1:N){
      # likelihood
      y[i] ~ dnorm(mu, 1/15^2) # normal likelihood with known SD of 15
    }
    # prior
    mu ~ dnorm(prior_m, 1/prior_sd^2) # normal prior on mu

  prior_m <- 100
  prior_sd <- 10
}
"

# create the JAGS model. n.adapt 'tunes' the samplers and don't contribute to the mcmc chains - rjags will warn you if n.adapt is too small
jags <- jags.model(file = textConnection(model), data = data_list, n.chains = 4, n.adapt = 1000)

update(jags, 1000) # warm up

samples = coda.samples(model = jags, variable.names = "mu", n.iter = 1000)

summary(samples)

effectiveSize(samples)

# look at the chains and the distribution of mu
plot(samples)

# look at convergence (Rhat) and autocorrelation
gelman.diag(samples)

par(mfrow=c(2,2))
autocorr.plot(samples, auto.layout = F) # produces a plot for each chain
par(mfrow=c(1,1))


# suppose that we do not know the SD and want to estimate it as well

# we need a new model specification (data stays the same)

model2 = "
  model {
    for (i in 1:N){
      # likelihood
      y[i] ~ dnorm(mu, 1/sigma^2) # normal likelihood with known SD of 15
    }
    # prior
    mu ~ dnorm(100, 1/90^2) # normal prior on mu
    sigma ~ dgamma(1, .1) # a gamma prior on SD
}
"

jags2 <- jags.model(file = textConnection(model2), data = data_list, n.chains = 4, n.adapt = 1000)

update(jags2, 1000) # warm up

samples2 = coda.samples(model = jags2, variable.names = c("mu", "sigma"), n.iter = 1000)

plot(samples2)
gelman.diag(samples2)

par(mfcol=c(2,4))
autocorr.plot(samples2, auto.layout = F)
par(mfcol=c(1,1))

# plot the joint posterior samples
plot(as.matrix(samples2), col = 'grey', type='l')




