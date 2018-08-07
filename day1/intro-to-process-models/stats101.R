####################### Probability Density ##################################

x <- seq(-3, 3, .01)
y <- dnorm(x = x, mean = 0, sd = 1)

plot(x, y, type = "l")



x <- seq(-40, 20, .01)
y <- dnorm(x = x, mean = -10, sd = 10)

plot(x, y, type = "l")


samp <- rnorm(10000, -10, 10)
hist(samp)

qnorm(p = .5, 0, 1)
