
# maxmimum likelihood modeling of Zhang and Luck (2008, nature)

de = read.table("zhang-luck08.dat")
# this data was uploaded by van den berg (2014) psych rev (paper here: http://www.cns.nyu.edu/malab/static/files/publications/2014%20Van%20den%20Berg%20Ma%20Psych%20Review.pdf)
# data here: http://www.cns.nyu.edu/malab/static/files/dataAndCode/Van_den_Berg_2014_DATA.zip
# the Zhang & Luck experiment is under E2 (see paper for other data sets)

str(de)
# participants saw displays of 1,2,3,6 items (setsize)
# error = distance from target in radians
# we add an index to tell us which setsize was tested on a given trial (1:4)
de$setsize_i = as.numeric(as.factor(de$setsize))

# the package 'circular' gives us circular probability distributions (those with support [-pi, pi])
# the von mises is the circular equivalent of the normal

# the function below gives the density for a mixture of a vonmises and uniform
# pmem is the probability the response comes from memory
de_mixture <- function(y, sd, pmem, mu=0, log=T){
  # delayed estimation mixture
  dvon = suppressWarnings(dvonmises(circular(y), mu = mu, kappa = 1/sd^2)) # suppresses messages about converting data to the circular class 
  p <- pmem*dvon + (1-pmem)*(1/(2*pi))
  if (log){
    return(log(p))
  }else{
    return(p)
  }
}

# the line below plots the expected distribution for values of sd and pmem
# try testing out different values
curve(de_mixture(x, sd = .5, pmem = .8, log=F), from=-pi, to=pi, ylim=c(0,.8), ylab='', xlab="Error (radians)")


### FUNCTIONS FOR MLE
# the three functions below use de_mixture to implement 3 different models

# zhang and luck estimated a different pmem and sd for each set size
zhang_luck <- function(par, y, N_i){
  # length(par) = 8 (4*sd + 4*pmem)
  N_n = length(unique(N_i))
  ll = 0
  for (i in 1:length(y)){
    ll = ll + de_mixture(y = y[i], sd = par[N_i[i]], pmem = par[N_n + N_i[i]], mu=0, log=T)
  }
  return(ll)
}

## test
# par = runif(8, min = 0, max = rep(c(20, 1), each = 4))
# zhang_luck(par = par, y = de$error[de$ppt==1], N_i = de$setsize_i[de$ppt==1])

# this version estimates k (number of items in memory) directly to determine pmem for a given setsize
zhang_luck_k <- function(par, y, N, N_i){
  # note this function takes the extra argument N
  # length(par) = 5 (4*sd + k)
  N_n = length(unique(N_i))
  ll = 0
  for (i in 1:length(y)){
    pm = min(par[N_n + 1]/N[i], 1)
    ll = ll + de_mixture(y = y[i], sd = par[N_i[i]], pmem = pm, mu=0, log=T)
  }
  return(ll)
}

## test
# par = runif(5, min = 0, max = c(rep(20, 4), 4))
# zhang_luck_k(par = par, y = de$error[de$ppt==1], N = de$setsize[de$ppt==1], N_i = de$setsize_i[de$ppt==1])

# a SDT "resource" model (no guessing)
# zhang and luck considered a very basic SDT model
# far more elaborate resource models have been proposed since (see, e.g., van den Berg et al., 2014, psych rev)

resource <- function(par, y, N_i){
  # length(par) = 4 (4*sd)
  N_n = length(unique(N_i))
  ll = 0
  for (i in 1:length(y)){
    ll = ll + de_mixture(y = y[i], sd = par[N_i[i]], pmem = 1, mu=0, log=T)
  }
  return(ll)
}

## test
# par = runif(4, min = 0, max = rep(20, 4))
# resource(par = par, y = de$error[de$ppt==1], N_i = de$setsize_i[de$ppt==1])

##### TASK ----
### Use the functions above to fit the three models to each participant's data (N = 8)
### Which provides the best fit for each individual?

