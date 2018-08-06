
library(rjags)
library(coda)

cd = read.table(file = "day2/bayesian-models/jags-change-det/rouder08-longdata-0.5.dat")

head(cd)

# ischange = was the trial a change trial? 1 = yes, 0 = no
# respchange = the number of trials that the participant said 'change' to
# ntrials = number of trials of that type for that participant in that condition

# so
# for ischange = 1 all the respchange responses are hits (ntrials - respchange = misses)
# for ischange = 0 all the respchange responses are false-alarms (ntrials - respchange = correct rejections)

k_model = "
  model {
    for (i in 1:n){
      y[i] ~ dbin(y.hat[i], ntrials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # p(resp = change)
      P[i] <- ifelse(ischange[i] == 1, 
                     a[i]*(d[i]+(1-d[i])*g[i]) + (1-a[i])*g[i], # p(hit)
                     a[i]*(1-d[i])*g[i] + (1-a[i])*g[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # 'Mass-at-chance' transformation
      kappa[i] <- K_s[id[i]]
      logit(a[i]) <- A_s[id[i]] # logit transformation
      logit(g[i]) <- G_s[id[i]]
    }
    # sample participant parameters from population distributions
    # in this model the distributions are independent but we could 
    # model correlations between the model parameters if they were 
    # sampled from, e.g., a multivariate normal
    for (s in 1:S){
      K_s[s] ~ dnorm(K_mu, 1/K_sigma^2)
      A_s[s] ~ dnorm(A_mu, 1/A_sigma^2)
      G_s[s] ~ dnorm(G_mu, 1/G_sigma^2)
    }
    
    K_mu ~ dnorm(3, 1/4^2) # k is typically around 3 but the sd is broad on this scale
    A_mu ~ dnorm(2.2, 1/4^2) 
    # ^ this is on the log odd scale. 3 implies the expectation that 
    # participants will pay attention on ~ 90% of the trials
    # again the prior is broad on this scale
    G_mu ~ dnorm(0, 1/4^2) # 0 = 0.5 on the probability scale
    
    K_sigma ~ dgamma(shape, rate)
    A_sigma ~ dgamma(shape, rate)
    G_sigma ~ dgamma(shape, rate)
    
    shape <- 1.01005 # mode = .1, SD = 10 (v. vauge)
    rate <- 0.1005012
  }
"

cd_list = list(
  "y" = cd$respchange,
  "ischange" = cd$ischange,
  "ntrials" = cd$ntrials,
  "N" = cd$N,
  "n" = nrow(cd),
  "S" = length(unique(cd$ppt)),
  "id" = cd$ppt
)


# initialize the jags model
k_jags <- jags.model(file = textConnection(k_model), data = cd_list, n.chains = 4, n.adapt = 1000)

# warm up the chains
update(k_jags, 1000)

params = c("K_mu", "A_mu", "G_mu", "K_sigma", "A_sigma", "G_sigma")
k_samples = coda.samples(model = k_jags, variable.names = params, n.iter = 1000)

summary(k_samples)

# a version of the model with a different k parameter for each set size

# to fit this we need to add some extra stuff to the data list
cd_list$N_i = as.numeric(as.factor(cd_list$N)) # an index (1:3) for which set size were are working with on trial i
cd_list$N_n = length(unique(cd_list$N)) # the number of different set sizes (3)

vary_k_model = "
  model {
    for (i in 1:n){
      y[i] ~ dbin(y.hat[i], ntrials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # p(resp = change)
      P[i] <- ifelse(ischange[i] == 1, 
                     a[i]*(d[i]+(1-d[i])*g[i]) + (1-a[i])*g[i], # p(hit)
                     a[i]*(1-d[i])*g[i] + (1-a[i])*g[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # 'Mass-at-chance' transformation
      
      # in this model individual capacities are stored in a matrix
      # row = individual, column = set size
      kappa[i] <- K_s[id[i], N_i[i]]
      logit(a[i]) <- A_s[id[i]] # logit transformation
      logit(g[i]) <- G_s[id[i]]
    }
    
    for (s in 1:S){
      for (ss in 1:N_n){
        # sample individual Ks from a normal with shared variance
        # but different mean for each set size
        K_s[s, ss] ~ dnorm(K_mu[ss], 1/K_sigma^2)
      }
      A_s[s] ~ dnorm(A_mu, 1/A_sigma^2)
      G_s[s] ~ dnorm(G_mu, 1/G_sigma^2)
    }
    
    for (ss in 1:N_n){
      K_mu[ss] ~ dnorm(3, 1/4^2)
    }
    A_mu ~ dnorm(2.2, 1/4^2) 
    G_mu ~ dnorm(0, 1/4^2)
    
    K_sigma ~ dgamma(shape, rate)
    A_sigma ~ dgamma(shape, rate)
    G_sigma ~ dgamma(shape, rate)
    
    shape <- 1.01005 # mode = .1, SD = 10 (v. vauge)
    rate <- 0.1005012
  }
"

# initialize the jags model
vary_k_jags <- jags.model(file = textConnection(vary_k_model), data = cd_list, n.chains = 4, n.adapt = 1000)

# warm up the chains
update(vary_k_jags, 1000)

params = c("K_mu", "A_mu", "G_mu", "K_sigma", "A_sigma", "G_sigma")
vary_k_samples = coda.samples(model = vary_k_jags, variable.names = params, n.iter = 1000)

summary(vary_k_samples)

gelman.diag(vary_k_samples)

plot(vary_k_samples[,"K_mu[1]"])
# the mean for set size = 2 has not converged - why do you think this is?

# re-fitting with more samples helps...
vary_k_samples = coda.samples(model = vary_k_jags, variable.names = params, n.iter = 10000)
gelman.diag(vary_k_samples)

# but this model probably should be reparameterized. If you're interested in how this could be done, please ask one of us

### Comparing these models with DIC ----

# DIC is like AIC and BIC but for hierarchical models https://en.wikipedia.org/wiki/Deviance_information_criterion
# it essentially penalizes the model for the 'effective' number of parameters it has (how this is estimated is tricky. You can't just could the number of paramaters in a hierarchical model)

# we can get it via the dic.samples function
DIC_k_jags = dic.samples(model = k_jags, n.iter = 1000, type = "pD")
DIC_vary_k_jags = dic.samples(model = vary_k_jags, n.iter = 1000, type = "pD")

DIC_k_jags
DIC_vary_k_jags

diffdic(DIC_k_jags, DIC_vary_k_jags)

# DIC is smaller for the fixed k version

### Posterior predictive samples ----

# we need to extract the samples to a matrix and then create a function to 
# turn the parameters into new data
k_samples_mat = as.matrix(k_samples)

k_ppsamples = function(mat, N, trials = 30){
  # takes the posterior samples and produces posterior predictive samples for 
  # a particular set size
  
  # sample new ks from the population distribution
  kappa = rnorm(nrow(mat), mean = mat[,'K_mu'], sd = mat[,'K_sigma'])
  k = ifelse(kappa > 0, kappa, 0) # this is the same as k = max(kappa, 0)
  
  d = k/N
  d = ifelse(d > 1, 1, d) # this is the same as d = min(k/N, 1)
  
  # sample new a and gs from their population distributions
  logita = rnorm(nrow(mat), mean = mat[,'A_mu'], sd = mat[,'A_sigma'])
  a = plogis(logita) # transform this parameter to its natural scale [0,1]
  
  logitg = rnorm(nrow(mat), mean = mat[,'G_mu'], sd = mat[,'G_sigma'])
  g = plogis(logita)
  
  h = a*(d + (1-d)*g) + (1 - a)*g
  f = a*(1 - d)*g + (1 - a)*g
  
  # simulate new data from a binomial
  h_rep = rbinom(n = length(h), size = trials, prob = h)/trials
  f_rep = rbinom(n = length(f), size = trials, prob = f)/trials
  
  return(cbind(f = f_rep, h = h_rep))
}

# these are the posterior predictive false-alarm and hit rates for different set sizes
pp_samples_N2 = k_ppsamples(mat = k_samples_mat, N = 2)
pp_samples_N5 = k_ppsamples(mat = k_samples_mat, N = 5)
pp_samples_N8 = k_ppsamples(mat = k_samples_mat, N = 8)

# let's look at quantiles
apply(pp_samples_N2, 2, FUN = quantile, prob=c(.025, .975))
apply(pp_samples_N5, 2, FUN = quantile, prob=c(.025, .975))
apply(pp_samples_N8, 2, FUN = quantile, prob=c(.025, .975))

### HARD task ----

## try ploting the data and the posterior predictive samples

# observed hit and false alarm rates can be found by 
cd$rate = with(cd, respchange/ntrials)


