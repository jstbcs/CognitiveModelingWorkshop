
# MLE Rouder et al (2008) PNAS

# get the MLE functions from the group script
source("day1/change-detection/mle-rouder08-group.R")

# the data is also read in under cd

head(cd)



# examine fit statistics
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