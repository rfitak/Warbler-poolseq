# Finding outlier windows using _GenWin_
We used the R package [_GenWin_](https://cran.r-project.org/web/packages/GenWin) ([Beissinger et al. 2015](https://doi.org/10.1186/s12711-015-0105-9)) to generate the windows to use for the analysis.  Rather than selecting an arbitrary window size (e.g., 25 kb) and step, GenWin fits a smoothed spline function to a metric (e.g., F<sub>ST</sub>) on each chromosome and subsequently uses changes in sign of the second derivative, or inflection points, to define windows. This directly taken from [Beissinger et al. 2015](https://doi.org/10.1186/s12711-015-0105-9):
>"We introduce a method for defining windows based on statistically guided breakpoints in the data, as a foundation for the analysis of multiple adjacent data points. This method involves first fitting a cubic smoothing spline to the data and then identifying the inflection points of the fitted spline, which serve as the boundaries of adjacent windows. This technique does not require prior knowledge of linkage disequilibrium, and therefore can be applied to data collected from individual or pooled sequencing experiments. Moreover, in contrast to existing methods, an arbitrary choice of window size is not necessary, since these are determined empirically and allowed to vary along the genome."

When using GenWin, we noticed there is a bug/issue in the current version v1.0.  If the `splineAnalyze` function only detects one inflection point (i.e. two windows), it throws an error.  I fixed this issue in the source code and created a new version: [GenWin v1.1](./data/GenWin_1.1.tar.gz).

_R code to run GenWin on each scaffold and append to the output file `genwin.out.csv`_

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
