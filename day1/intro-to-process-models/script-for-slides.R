############### First High-Thrshold Model #################################################

#negative log likelihood of high-threshold model
nll.ht <- function(par, y){ #1. argument: vector with parameters, 2. arg: vector with data (h, m, f, c)
  d <- par[1]
  g <- par[2]
  p <- 1:4 # reserve space
  p[1] <- d + (1 - d) * g   #probability of a hit
  p[2] <- 1 - p[1]          # probability of a miss
  p[3] <- g                 # probability of a false alarm
  p[4] <- 1 - p[3]          #probability of a correct rejection
  return(-sum(y * log(p)))
}

y <- c(75, 25, 30, 20) #h, m, f, c
par <- c(.5, .5) #starting values for probability parameters
out <- optim(par, nll.ht, y = y)
out$par
out$value

################ Model Comparison High-Threshold Model ####################################


#negative log likelihood for any one condition
nll.condition <- function(par, y){ #assign par=c(d, g), y = c(h, m, f, c)
  p <- 1:4
  d <- par[1]
  g <- par[2]
  p[1] <- d + (1 - d) * g
  p[2] <- 1 - p[1]
  p[3] <- g
  p[4] <- 1 - p[3]
  return(-sum(y * log(p)))
}

#negative log likelihood for General Model:
#assign par4 = (d1, g1, d2, g2), y8 = (h1, m1, f1, c1, h2, m2, f2, c2)
nll.g <- function(par4, y8){            
  nll.condition(par4[1:2], y8[1:4]) +     #condition 1
    nll.condition(par4[3:4], y8[5:8])       #condition 2
}

#negative log likelihood for Model 1:
#assign par3 = (d, g1, g2), y8 = (h1, m1, f1, c1, h2, m2, f2, c2)
nll.1 <- function(par3, y8){            
  nll.condition(par3[1:2], y8[1:4]) +     #condition 1
    nll.condition(par3[c(1, 3)], y8[5:8])       #condition 2
}

#negative log likelihood for Model 2:
#assign par3 = (d1, d2, g), y8 = (h1, m1, f1, c1, h2, m2, f2, c2)
nll.2 <- function(par3, y8){            
  nll.condition(par3[c(1, 3)], y8[1:4]) +     #condition 1
    nll.condition(par3[2:3], y8[5:8])       #condition 2
}



### Data analysis
dat <- c(22, 28, 22, 28   #h, m, f, c for condition 1
         , 35, 15, 21, 29) #h, m, f, c for condition 2

#General Model
par.m <- c(.5, .5, .5, .5) #starting values
mod.g  <- optim(par.m, nll.g,y8 = dat, hessian = T)

#Model 1
par.m <- c(.5, .5, .5) #starting values
mod.1 <- optim(par.m, nll.1, y8 = dat, hessian = T)

#Model 2
par.m <- c(.5, .5, .5) #starting values
mod.2 <- optim(par.m, nll.2, y8 = dat, hessian = T)



### Model comparison

G1 <- 2*(mod.1$value - mod.g$value)
G2 <- 2*(mod.2$value - mod.g$value)

qchisq(.95, df = 1) #Critical value for alpha = .05

c(m1 = 1 - pchisq(G1, df = 1), m2 = 1 - pchisq(G2, df = 1)) #p-values