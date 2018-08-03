
library(plyr)
library(RCurl)

intext=getURL("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/wmPNAS2008/lk2data.csv")
data=read.csv(text=intext)

head(data)

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

