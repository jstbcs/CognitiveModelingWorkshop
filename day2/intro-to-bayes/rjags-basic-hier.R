
library(rjags)
library(coda)

# a basic hierarchical model where we estimate a mean for each participant
# as coming from a population distribution of means

#### SIMULATE DATA -----

I = 30 # number of participants
J = 10 # number of trials per participant
N = I*J # number of observations

set.seed(123)
# sample the 'true' participant means from a normal
# with a population mean of 700 and sd of 100
mu_i = rnorm(I, 700, 100)

# create the data frame to hold the simulated data
hier_data = expand.grid(ppt = 1:I, trial = 1:J, y = NA)
for (i in 1:I){
  y_i = rnorm(J, mu_i[i], 100) # generate the trials for this participant
  hier_data$y[hier_data$ppt == i] <- y_i # and add to the data frame
}

# observed participant means
(xbars = aggregate(y ~ ppt, data = hier_data, FUN = mean))

#### MODEL -----

# JAGS model code for the basic hierarchical model
model = "
  model {
    for (n in 1:N){
      # likelihood
      y[n] ~ dnorm(mu[id[n]], 1/sigma^2)
    }
    # sample participant parameters from population distribution
    for (i in 1:I){
      mu[i] ~ dnorm(m, 1/s^2)
    }
    # priors
    m ~ dnorm(800, 1/200^2) 
    s ~ dgamma(2, .01)
    sigma ~ dgamma(2, .01) # a gamma prior on SD
}
"

head(hier_data)

# put the data into a list
data_list = list(
  'y' = hier_data$y,
  'id' = hier_data$ppt, # participant id for each observation
  'I' = length(unique(hier_data$ppt)), # number of participantd
  'N' = nrow(hier_data) # number of observations
)

# initialize the jags model
jags <- jags.model(file = textConnection(model), data = data_list, n.chains = 4, n.adapt = 1000)

# warm up the chains
update(jags, 1000)

samples = coda.samples(model = jags, variable.names = c("mu", "sigma", "m", "s"), n.iter = 1000)

plot(samples[,"m"])
plot(samples[,"s"])
plot(samples[,"sigma"])

gelman.diag(samples)

# plot histogram of samples for some participants
samples_matrix = as.matrix(samples) # convert mcmc.list into a matrix

head(samples_matrix)

plot_ids = sample(1:I, 8) # select 8 random participant numbers

par(mfrow=c(4,2), mar=c(4,4,1,1)) # set layout for plots
for (i in plot_ids){
  hist(samples_matrix[,paste0("mu[", i, ']')], main = paste0("ppt = ", i), xlab="", col="lightgreen")
}

# calculate crucial quantities of interest
# posterior means (apply takes a matrix and applies a function to the rows [MARGIN=1] or columns [MARGIN=2])
(post_m <- apply(samples_matrix, MARGIN = 2, FUN = mean))
# posterior medians
apply(samples_matrix, MARGIN = 2, FUN = median)
# 95% credible intervals
(post_cis <- apply(samples_matrix, MARGIN = 2, FUN = quantile, probs = c(.025, .975)))

# you can also get this information via
summary(samples)

# lets plot the posterior means against the 'true' means
post_m = post_m[2:31]
post_cis = post_cis[,2:31]

par(mfrow=c(1,1), mar=c(5,4,3,2)) # back to a 1 panel plot

plot(x = mu_i, y = post_m, xlab="True", ylab="Estimated (95% CI)", main="True means vs estimated", pch=16, col="red")
segments(x0 = mu_i, y0 = post_cis[1,], x1 = mu_i, y1 = post_cis[2,]) # error bars
abline(0,1, lty=2, col="grey") # line of 1:1 correspondence

# plot the posteriors against the sample means
plot(x = xbars$y, y = post_m, xlab=bquote(bar(y)[i]), ylab="Estimated (95% CI)", main="Participant means vs estimated", pch=16, col="red")
segments(x0 = xbars$y, y0 = post_cis[1,], x1 = xbars$y, y1 = post_cis[2,]) # error bars
abline(0,1, lty=2, col="grey") # line of 1:1 correspondence

## Posterior predictive check

hier_ppsamples <- function(m, s, sigma){
  # this function takes the hyperparameters from the above model
  # and generates new observations
  m_n = rnorm(length(m), m, s) # samples means from population distribution
  
  y_rep = rnorm(length(m), m_n, sigma) # generate new "data"
  
  return(y_rep)
}

pp_samples = hier_ppsamples(samples_matrix[,'m'], samples_matrix[,'s'], samples_matrix[,'sigma'])

par(mfrow=c(1,2))
hist(pp_samples, main="Posterior Predictive Samples", probability = T, xlab="y", breaks=30, col="lightblue")

hist(hier_data$y, xlim=range(pp_samples), main="Data", xlab='y', col="lightgreen")
