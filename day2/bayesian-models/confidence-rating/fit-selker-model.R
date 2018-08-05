
library(rjags)
library(coda)
# this loads the data from Pratte et al (2010) Separating Mnemonic Process From Participant and Item Effects in the Assessment of ROC Asymmetries. JEPLMC, 36(1), 224.
# reshaped with code provided by Selker et al. https://osf.io/b8gcf/

load("day2/bayesian-models/confidence-rating/pratte10.RData")

str(pratte10_list)

# xNoise = data from new trials
# xSignal = data from old trials
# xNoise and xSignal are 97x240 matricies. Rows = participants, columns = trials 
# observations are ratings 1:6. 
# 1 = sure new, 2 = believe new, 3 = guess new, 4 = guess old, 5 = believe old, 6 = sure old
# nNoise, nSignal = number of noise and signal trials, respectively, per participant
# nCat = number of confidence categories (6)
# nSubjs = N participants (97)

# the model ("HierSDT_model.txt") comes from:
# Selker et al. Parsimonious Estimation of Signal Detection Models from Confidence Ratings. Pre-print available here: https://osf.io/b6z8e/

## It can take a long time to fit the model to all participants so instead we can
## select a subset of participants

SUBSET = TRUE # set to false if you want to fit to all participants

if (SUBSET){
  nselect = 20 # how many do we want
  ppts = sample(1:97, size = nselect)
  
  pratte10_list$xNoise = pratte10_list$xNoise[ppts,]
  pratte10_list$xSignal = pratte10_list$xSignal[ppts,]
  
  pratte10_list$nNoise = pratte10_list$nNoise[ppts]
  pratte10_list$nSignal = pratte10_list$nSignal[ppts]
  
  pratte10_list$nSubjs = nselect
}


# initialize the jags model
rating_jags <- jags.model(file = "day2/bayesian-models/confidence-rating/HierSDT_model.txt", data = pratte10_list, n.chains = 4, n.adapt = 1000)

# warm up the chains
update(rating_jags, 1000)

params <- c("muMu", "mu", "sigmaMu", "sigma", "aMu", "a", "bMu", "b")
rating_samples = coda.samples(model = rating_jags, variable.names = params, n.iter = 1000)

gelman.diag(rating_samples)

plot(rating_samples)
