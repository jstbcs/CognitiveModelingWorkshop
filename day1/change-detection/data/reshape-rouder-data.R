
library(plyr)
library(RCurl)

intext=getURL("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/wmPNAS2008/lk2data.csv")
data=read.csv(text=intext)

head(data)

## Wide format for MLE

counts = ddply(data, c('sub', 'prch', 'N'), summarize,
      H = sum(ischange == 1 & resp == 1),
      M = sum(ischange == 1 & resp == 0),
      Fa = sum(ischange == 0 & resp == 1),
      Cr = sum(ischange == 0 & resp == 0))

counts$N = paste0("N", counts$N)

counts_wide = 
counts %>% 
  gather(variable, value, -(sub:N)) %>%
  unite(temp, N, prch, variable) %>%
  spread(temp, value)

colorder = c()
for (i in c(0.3, 0.5, 0.7)){
  for (j in c("N2", "N5", "N8")){
    colorder <- c(colorder, paste(j, i, c("H", "M", "Fa", "Cr"), sep="_"))
  }
}

# re-order columns
counts_wide = counts_wide[, colorder]
apply(counts_wide, 1, sum)

write.table(x = counts_wide, file = "rouder08-data-full.dat")

# only the 50:50 trials
counts_wide_0.5 = counts_wide[,grep(colorder, pattern = "0.5")]

write.table(x = counts_wide_0.5, file = "rouder08-data-0.5.dat")

# -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ -~ 
## Long format for JAGS

data_long = ddply(data, c("sub", "prch", "N", "ischange"), summarize,
      respchange = sum(resp), ntrials = length(resp))

colnames(data_long)[1] = "ppt"

data_long$ppt = as.numeric(as.factor(data_long$ppt)) # renumber participants 1:23

setwd("../../../day2/bayesian-models/jags-change-det/")

write.table(x = data_long, file = "rouder08-longdata-full.dat")

data_long_0.5 = subset(data_long, prch==0.5)

write.table(x = data_long_0.5, file = "rouder08-longdata-0.5.dat")
