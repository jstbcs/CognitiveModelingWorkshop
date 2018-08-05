####################### Probability Density ##################################

x <- seq(-3, 3, .01)
y <- dnorm(x = x, mean = 0, sd = 1)

plot(x, y, type = "l")