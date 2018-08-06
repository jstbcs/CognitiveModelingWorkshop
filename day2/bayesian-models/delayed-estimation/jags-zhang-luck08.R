
library(rjags)
library(coda)

# we need the von mises distribution for this example
# this doesn't come as standard but can be added to JAGS
# for instruction on how to get this, see:
# see https://github.com/yeagle/jags-vonmises and https://link.springer.com/article/10.3758%2Fs13428-013-0369-3
load.module("vonmises")
load.module("dic")

# Hierarchical models for Zhang and Luck (2008)

de = read.table("day2/bayesian-models/delayed-estimation/zhang-luck08.dat")

head(de)

# the model fit by zhang and luck estimated a separate P_m and sd for each set size
# and each participant!
# the code for the hierarchical version can be found here: "day2/bayesian-models/delayed-estimation/models/z-l-mixture.txt"
# it runs very slowly (if at all) and probably needs to be reparameterized

# below we fit a model where a capacity (k) is estimated for each participant
# as well as a sd for each set size

# first we fit the 'resource' model which estimates a different 
# sd for each set size

de_list = list(
  mu = pi, # the vonmises distribution supports [0, 2*pi] 
  y = de$error + pi, # so we have to add pi to the data to re-center
  N_i = as.numeric(as.factor(de$setsize)), # index of set size condition
  N_n = length(unique(de$setsize)), # number of set sizes
#  N = de$setsize,
  id = de$ppt,
  S = length(unique(de$ppt)),
  n = nrow(de)
)

# initialize the jags model
resource_jags <- jags.model(file = "day2/bayesian-models/delayed-estimation/models/z-l-resource.txt", data = de_list, n.chains = 2, n.adapt = 1000)

# warm up the chains
update(resource_jags, 1000)

params = c("SD_s", "SD_mu", "SD_sigma")
resource_samples = coda.samples(model = resource_jags, variable.names = params, n.iter = 1000)

gelman.diag(resource_samples)

plot(resource_samples[,'SD_mu[1]'])
plot(resource_samples[,'SD_mu[2]'])
plot(resource_samples[,'SD_mu[3]'])
plot(resource_samples[,'SD_mu[4]'])

### Task
# the sd samples are exp transformed
# could you transform them back to their natural scale and plot the histograms?

## Mixture model with k limit

de_list$N = de$setsize # add set size to data list

# initialize the jags model
mixture_k_jags <- jags.model(file = "day2/bayesian-models/delayed-estimation/models/mixture_k.txt", data = de_list, n.chains = 2, n.adapt = 1000)

# warm up the chains
update(mixture_k_jags, 1000)

params = c("K_s", "K_mu", "K_sigma", "SD_s", "SD_mu", "SD_sigma")
mixture_k_samples = coda.samples(model = mixture_k_jags, variable.names = params, n.iter = 1000)

plot(mixture_k_samples[,'K_mu'])
plot(mixture_k_samples[,'SD_mu[1]'])
plot(mixture_k_samples[,'SD_mu[2]'])
plot(mixture_k_samples[,'SD_mu[3]'])
plot(mixture_k_samples[,'SD_mu[4]'])
