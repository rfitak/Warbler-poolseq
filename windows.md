# Finding outlier windows using _GenWin_
xxx

```R
# Load Library (v1.1)
library(GenWin)

# Read in grenedalf Fst
data <- read.csv("grenedalf.fst.csv", header = T)
   #14,491,335 sites

# Trim to just B1 vs B2 Fst without nan
data2 <- subset(data, !is.na(B1.clean.sorted.B2.clean.sorted.fst))

# Get full list of scaffolds
scaffolds <- unique(data$chrom)
   #29,163 scaffolds in list

for (i in 1:length(scaffolds)){
   tmp <- subset(data2, chrom  == scaffolds[i])

   # This "if" statement avoids an error int he # of SNPs
   # The # SNPs must be >2*N+1 (N = 2, or the order of the polynomial in smooth.pspline)

   if (nrow(tmp) <= 5){
      message(paste0("Skipped scaffold ", i, " which had only ", nrow(tmp), " SNPs."))
   } else {
   win <- splineAnalyze(Y = tmp$B1.clean.sorted.B2.clean.sorted.fst, 
   map = tmp$start,
   smoothness = 100,
   plotRaw = FALSE,
   plotWindows = FALSE,
   method = 4)

   # Add the scaffold name to each window under the column "Chrom"
   win <- cbind(Chrom = rep(scaffolds[i], nrow(win$windowData)), win$windowData)

   # Write to output table, append each scaffold but don't reprint the header each time.
   write.table(win,
      file = "genwin.out.csv",
      append = T,
      quote = F,
      sep = ",",
      row.names = F,
      col.names = !file.exists("genwin.out.csv"))
   message(paste0("Finished scaffold ", i))
   }
}
```
