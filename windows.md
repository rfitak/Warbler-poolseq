# Finding outlier windows using _GenWin_
xxx

```R
# R code
library(GenWin)

# Read in grenedalf Fst
data <- read.csv("grenedalf.fst.csv", header = T)
   #14,491,335 sites

# Trim to just B1 vs B2 Fst without nan
data2 <- subset(data, !is.na(B1.clean.sorted.B2.clean.sorted.fst))

# Get full list of scaffolds
scaffolds <- unique(data$chrom)
   #29,163 scaffolds in list

tmp <- subset(data2, chrom  == scaffolds[1])
win <- splineAnalyze(Y = tmp$B1.clean.sorted.B2.clean.sorted.fst, 
   map = tmp$start,
   smoothness = 100,
   plotRaw = FALSE,
   plotWindows = FALSE,
   method = 4)


 write.table(win$windowData, file = "delete.R.csv", append = T, quote = F, sep = ",", row.names = F, col.names = T)
```
