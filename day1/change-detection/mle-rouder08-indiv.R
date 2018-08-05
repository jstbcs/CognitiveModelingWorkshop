
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

sdt_fit$AIC
k_fit$AIC

sdt_fit$BIC
k_fit$BIC

#### FIT TO INDIVIDUALS ----

S = nrow(cd) # number of participants

# create matrices to hold the resulting parameter estimates
# 1 row per participant, 1 column per parameter
estimates_fix_k <- matrix(NA, nrow = S, ncol = 3)
colnames(estimates_fix_k) <- c("k", "a", "g")

estimates_vary_k <- matrix(NA, nrow = S, ncol = 5)
colnames(estimates_fix_k) <- c("k1", "k2", "k3", "a", "g")

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
