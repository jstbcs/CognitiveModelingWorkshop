## Get help

help.start()
citation()
sessionInfo()

?sum
?mean
help("sd")

??correlation

## Basic calculator

1+1
12^4
log(23)

x <- 10
y <- 6
x + y

z <- x + y
y/z

#Warning: R is case sensitive
X <- 5 # Does not replace variable x
x

## Variables

name <- "Jed Bartlet"  # Character
age <- 65           # Integer / Double
president <- TRUE         # Boolean

staff <- c("Josh", "CJ", "Sam")
staff

staff[2]
staff[c(1, 3)]
staff <- staff[c(1, 3)]
staff[2:3]

#Matrix
m <- matrix(c(1:9), nrow = 3)
m[1, ] # First row
m[, 3] # Third column
m[1, 3] # Element in the first row and third column

#Dataframe
name <- c("Ben", "Kathy", "Paul", "Julia", "Jeff")
age <- c(24, 43, 32, 27, 60)
job <- rep("teacher", 5)
friends <- c(5, 2, 0, 3, 6)
pets <- c(0, 2, 4, 1, 1)
income <- c(":/", ":)", ":(", ":(", ":)")

teachers <- data.frame(name, age, job, friends, pets, income)
save(teachers, file = "day1/intro-to-R/data/teachers.RData")

teachers$ratio <- teachers$pets / teachers$friends
teachers


## Using functions

mean(teachers$age[teachers$income == ":)"])
tapply(teachers$age, teachers$income, mean)

set.seed(666)
x <- rnorm(50, 10, 2)
y <- rnorm(50, 12, 4)
dat.fun <- data.frame(sub = 1:50, x, y)

write.csv(dat.fun, "day1/intro-to-R/data/example_fun.csv", row.names = F)

example.dat <- read.csv(file = "day1/intro-to-R/data/example_fun.csv")

head(dat.fun)

t.test(example.dat$x, example.dat$y, var.equal = T)

## Get data from the internet

daturl <- curl::curl("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/wmPNAS2008/lk2clean.csv")
dat <- read.csv(daturl, header = T)

head(dat)

#and clean it up
dat.pretty <- dat[, c("sub", "blk", "trl", "prch", "N", "ischange", "resp")]
dat.pretty$accuracy <- dat.pretty$ischange == dat.pretty$resp

dat.pretty$ischange <- factor(dat.pretty$ischange, labels = c("unchanged", "changed"))

head(dat.pretty)

mean.acc <- with(dat.pretty, tapply(accuracy, list(prch, N, ischange), mean))
mean.acc

#plot it
layout(matrix(1:2, ncol = 2))
matplot(mean.acc[,,1], main = "Unchanged", ylab = "Accuracy"
        , type = "b", pch = colnames(mean.acc), xaxt = "n")
axis(1, at = 1:3, labels = rownames(mean.acc))
matplot(mean.acc[,,2], main = "Changed", ylab = "Accuracy"
        , type = "b", pch = colnames(mean.acc), xaxt = "n")
axis(1, at = 1:3, labels = rownames(mean.acc))

matplot(mean.acc[,,1] # data
        , main = "Unchanged", ylab = "Accuracy" #labels: title and y-axis
        , type = "b", pch = colnames(mean.acc)) #type = "both" line and point, point type = set size

#count it
(tab <- table(dat.pretty$resp, dat.pretty$ischange))

(outcomes <- c("hits" = tab[2, 2], "misses" = tab[1, 2]
               , "fa" = tab[2, 1], "cr" = tab[1, 1]))

barplot(outcomes, col = "mediumvioletred")


## Write your own function

hello <- function(){
  return("Hello World!")
}

hello()

gimme.sum <- function(x, y){
  sum.xy <- x + y
  return(sum.xy)
}

gimme.sum(99, 567)

## New function
new.f <- function(x, y){
  sdy <- sd(y)
  val <- x/sdy
  return(val)
}

new.f(1:3, c(1, 1, 1))
