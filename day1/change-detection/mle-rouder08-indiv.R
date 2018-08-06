
# MLE Rouder et al (2008) PNAS

# get the MLE functions from the group script
source("day1/change-detection/mle-rouder08-group.R")

# the data is also read in under cd
head(cd)

# function to calculate fit statistics from -LL
fit_stats <- function(nLL, n, p){
  # nLL = negative log liklihood
  # n = number of observations
  # p = number of parameters
  
  deviance = 2*nLL
  aic = deviance + 2*p
  bic = deviance + p*log(n)
  
  return(list("D" = deviance, "AIC" = aic, "BIC" = bic))
}

sdt_fit = fit_stats(nLL = sdt_res$value, n = sum(group_data), p = 4)

k_fit = fit_stats(nLL = k_res$value, n = sum(group_data), p = 3)

vary_k_fit = fit_stats(nLL = vary_k_res$value, n = sum(group_data), p = 5)

sdt_fit$AIC
k_fit$AIC
vary_k_fit$AIC

sdt_fit$BIC
k_fit$BIC
vary_k_fit$BIC

#### FIT TO INDIVIDUALS ----

S = nrow(cd) # number of participants

# create matrices to hold the resulting parameter estimates
# 1 row per participant, 1 column per parameter
estimates_fix_k <- matrix(NA, nrow = S, ncol = 3)
colnames(estimates_fix_k) <- c("k", "a", "g")

estimates_vary_k <- matrix(NA, nrow = S, ncol = 5)
colnames(estimates_vary_k) <- c("k1", "k2", "k3", "a", "g")

estimates_sdt <- matrix(NA, nrow = S, ncol = 4)
colnames(estimates_sdt) <- c("d1", "d2", "d3", "c")

# create a matrix to hold the -log likelihood for each individual (row)
# and each model (col)
fit_statistics <- matrix(NA, nrow = S, ncol = 5)
colnames(fit_statistics) <- c("LL_vac", "LL_fix_k", "LL_vary_k", "LL_sdt", "N_obs")

# this loop takes the data from each row (participant) and fits the three models
for (s in 1:S){
  # get the data for this subject
  tmp.dat = as.integer(cd[s,])

  # model that freely estimates response frequencies
  fit_statistics[s,1] <- ll.vacuous(y = tmp.dat)
  
  # fixed k
  par = runif(n = 3, min = 0, max = c(max(N), 1, 1))
  k_res_s = optim(par, ll.fixed_k, y = tmp.dat)
  
  fit_statistics[s,2] <- k_res_s$value # add estimates and LL to matrices
  estimates_fix_k[s,] <- k_res_s$par
  
  # variable k
  par = runif(n = 5, min = 0, max = c(rep(max(N),3), 1, 1))
  vary_k_res_s = optim(par, ll.vary_k, y = tmp.dat)
  
  fit_statistics[s,3] <- vary_k_res_s$value
  estimates_vary_k[s,] <- vary_k_res_s$par
  
  ## sdt model
  par = runif(n = 4, min = 0, max = c(5, 5, 5, 5))
  sdt_res_s = optim(par, ll.sdt.ev, y = tmp.dat)
  
  fit_statistics[s,4] <- sdt_res_s$value
  estimates_sdt[s,] <- sdt_res_s$par
  
  fit_statistics[s,5] = sum(tmp.dat)
}
# remove stuff we no longer need...
rm(list = c("tmp.dat", "k_res_s", "vary_k_res_s", "sdt_res_s"))

# look at resulting parameter estimates
hist(estimates_fix_k[,'k'], main="Fixed k", xlab="k estimate")


#################### Model Comparison #######################

##Let's do AIC first
AIC.ind <- fit_statistics
for(s in 1:S){
  for(m in 1:M){
    AIC.ind[s, m] <- fit_stats(nLL = fit_statistics[s, m], n = fit_statistics[s, 5], p = npar[m])$AIC
  }
  AIC.ind[s, 5] <- order(AIC.ind[s, 1:4])[1]
}

colnames(AIC.ind) <- c("vac", "fix_k", "vary_k", "sdt", "winner")
AIC.ind

##BIC
BIC.ind <- fit_statistics
M <- ncol(BIC.ind)
npar <- c(12, 3, 5, 4)

for(s in 1:S){
  for(m in 1:M){
    BIC.ind[s, m] <- fit_stats(nLL = fit_statistics[s, m], n = fit_statistics[s, 5], p = npar[m])$BIC
  }
  BIC.ind[s, 5] <- order(BIC.ind[s, 1:4])[1]
}

colnames(BIC.ind) <- c("vac", "fix_k", "vary_k", "sdt", "winner")
BIC.ind


##################### More Stuff #####################################

#### Unequal Variance Signal Detection Model

## fit sdt model
par = runif(n = 7, min = 0, max = 3)
sdt_res_uv = optim(par, ll.sdt.uv, y = group_data)
sdt_res_uv$par

## fit sdt model
par = runif(n = 4, min = 0, max = 3)
sdt_res = optim(par, ll.sdt.ev, y = group_data)
sdt_res$par

c(sdt_res_uv$value, sdt_res$value)

## Try with differen random seeds set.seed(123)



##### Dealing with zero counts


# create a matrix to hold the -log likelihood for each individual (row)
# and each model (col)
fit_statistics <- matrix(NA, nrow = S, ncol = 5)
colnames(fit_statistics) <- c("LL_vac", "LL_fix_k", "LL_vary_k", "LL_sdt", "N_obs")

# this loop takes the data from each row (participant) and fits the three models
for (s in 1:S){
  # get the data for this subject
  tmp.dat = as.integer(cd[s,]) + .5
  
  # model that freely estimates response frequencies
  fit_statistics[s,1] <- ll.vacuous(y = tmp.dat)
  
  # fixed k
  par = runif(n = 3, min = 0, max = c(max(N), 1, 1))
  k_res_s = optim(par, ll.fixed_k, y = tmp.dat)
  
  fit_statistics[s,2] <- k_res_s$value # add estimates and LL to matrices
  estimates_fix_k[s,] <- k_res_s$par
  
  # variable k
  par = runif(n = 5, min = 0, max = c(rep(max(N),3), 1, 1))
  vary_k_res_s = optim(par, ll.vary_k, y = tmp.dat)
  
  fit_statistics[s,3] <- vary_k_res_s$value
  estimates_vary_k[s,] <- vary_k_res_s$par
  
  ## sdt model
  par = runif(n = 4, min = 0, max = c(5, 5, 5, 5))
  sdt_res_s = optim(par, ll.sdt.ev, y = tmp.dat)
  
  fit_statistics[s,4] <- sdt_res_s$value
  estimates_sdt[s,] <- sdt_res_s$par
  
  fit_statistics[s,5] = sum(tmp.dat)
}
# remove stuff we no longer need...
rm(list = c("tmp.dat", "k_res_s", "vary_k_res_s", "sdt_res_s"))

# look at resulting parameter estimates
hist(estimates_fix_k[,'k'], main="Fixed k", xlab="k estimate")

