
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
