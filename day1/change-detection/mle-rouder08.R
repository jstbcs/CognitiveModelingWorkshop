
# MLE Rouder et al (2008) PNAS

cd = read.table(file = "data/rouder08-data-0.5.dat")

group_data = apply(cd, 2, sum)

# the data frame gives numbers of hits, misses, false-alarms, and correct rejections
# for three set sizes: N = 2,5,8

N = c(2,5,8)
N_i = rep(1:length(N), each=4) # index

#Multinomial Negative Log-Likelihood
negLL <- function(y,p){
  a=ifelse(y==0 & p==0,0, y*log(p))
  -sum(a)
}

cowan_k <- function(k, a, g, N){
  d = min(1,k/N) # p(probe in memory)
  
  p = 1:4
  p[1] = a*(d+(1-d)*g)+(1-a)*g # p(hit)
  p[2] = 1-p[1] # p(miss)
  p[3] = a*(1-d)*g+(1-a)*g # p(false-alarm)
  p[4] = 1-p[3] # p(correct rejection)
  return(p)
}

sdt <- function(d, c, s){
  # this is a simplified version of the sdt 
  # model used by rouder et al.
  p = 1:4
  p[1] = pnorm((d - c)/s) # p(hit)
  p[2] = 1 - p[1] # p(miss)
  p[3] = pnorm(- c) # p(false-alarm)
  p[4] = 1 - p[3] # p(correct rejection)
  return(p)
}

# test this function out... plot an ROC curve
# m = sapply(seq(0,5, .1), FUN = function(x) sdt(d=0, c = x, s = 1))
# plot(m[3,], m[1,], type='l')

# Likelihood functions

ll.vacuous <- function(y){
  ll=0
  lenY = length(y)
  y1=y[(1:lenY)%%2==1]
  y2=y[(1:lenY)%%2==0]
  n=(rep((y1+y2),each=2))
  p=y/n
  ll = negLL(y, p)
  return(ll)
}

ll.fixed_k <- function(par, y){
  # length(par) == 3 (k, a, g)
  ll=0
  for(i in 1:length(N)){ # for each set size
    p = cowan_k(k = par[1], a = par[2], g = par[3], N = N[i])
    ll = ll + negLL(y[N_i==i], p)
  }
  if(any(c(par < rep(0,3), par > c(max(N),1,1)))){
    ll = ll + 10000 # penalty for going out of range
  }
  return(ll)
}

ll.sdt.ev <- function(par, y){
  # length(par) == 4 (d1, d2, d3, c)
  ll=0
  for(i in 1:length(N)){ # for each set size
    p = sdt(d = par[i], c = par[length(N)+1], s = 1)
    ll = ll + negLL(y[N_i==i], p)
  }
  return(ll)
}

# get LL from vacuous model 
ll.vacuous(y = group_data)

## fit k model
# starting values
par = runif(n = 3, min = 0, max = c(max(N), 1, 1))
k_res = optim(par, ll.fixed_k, y = group_data)
k_res$value

k_res$par

## fit sdt model
par = runif(n = 4, min = 0, max = c(5, 5, 5, 5))
sdt_res = optim(par, ll.sdt.ev, y = group_data)
sdt_res$value

sdt_res$par

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



