# Finding outlier windows using _GenWin_
We used the R package [_GenWin_](https://cran.r-project.org/web/packages/GenWin) ([Beissinger et al. 2015](https://doi.org/10.1186/s12711-015-0105-9)) to generate the windows to use for the analysis.  Rather than selecting an arbitrary window size (e.g., 25 kb) and step, GenWin fits a smoothed spline function to a metric (e.g., F<sub>ST</sub>) on each chromosome and subsequently uses changes in sign of the second derivative, or inflection points, to define windows. This directly taken from [Beissinger et al. 2015](https://doi.org/10.1186/s12711-015-0105-9):
>"We introduce a method for defining windows based on statistically guided breakpoints in the data, as a foundation for the analysis of multiple adjacent data points. This method involves first fitting a cubic smoothing spline to the data and then identifying the inflection points of the fitted spline, which serve as the boundaries of adjacent windows. This technique does not require prior knowledge of linkage disequilibrium, and therefore can be applied to data collected from individual or pooled sequencing experiments. Moreover, in contrast to existing methods, an arbitrary choice of window size is not necessary, since these are determined empirically and allowed to vary along the genome."

When using GenWin, we noticed there is a bug/issue in the current version v1.0.  If the `splineAnalyze` function only detects one inflection point (i.e. two windows), it throws an error.  I fixed this issue in the source code and created a new version: [GenWin v1.1](./data/GenWin_1.1.tar.gz).

### Step 1: Generate the windows
_R code to run GenWin on each scaffold and append to the output file `genwin.out.csv`_

```R
# Load Library (v1.1)
library(GenWin)

# Read in PoPoolation Fst file
data <- read.table("pools.fst", header = F, sep = "\t")[,c(1,2,6)]
colnames(data) <- c("Chrom", "Pos", "Fst")
   #14,439,829 sites

# Trim to just B1 vs B2 Fst without nan
data$Fst <- as.numeric(gsub("1:2=", "", data$Fst))
#data2 <- subset(data, !is.na(B1.clean.sorted.B2.clean.sorted.fst))

# Get full list of scaffolds
scaffolds <- unique(data$Chrom)
   #28,878 scaffolds in list

# Run GenWin by scaffold
for (i in 1:length(scaffolds)){
   tmp <- subset(data, Chrom == scaffolds[i])

   # This "if" statement avoids an error in the # of SNPs
   # The # SNPs must be >2*N+1 (N = 2, or the order of the polynomial in smooth.pspline)
   if (nrow(tmp) <= 5){
      message(paste0("Skipped scaffold ", scaffolds[i], " which had only ", nrow(tmp), " SNPs."))
   } else {
	   if (diff(range(tmp$Pos)) >= 100){
		   win <- splineAnalyze(Y = tmp$Fst,
                          map = tmp$Pos,
                          smoothness = 100,
                          plotRaw = FALSE,
                          plotWindows = FALSE,
                          method = 4)
	   } else if (diff(range(tmp$Pos)) < 100) {
              # This statement is to resolve an error when the SNPs on the scaffold
              # cover a region less than the smoothness parameter (100 bp)
                   win <- splineAnalyze(Y = tmp$Fst,
			  map = tmp$Pos,
			  smoothness = diff(range(tmp$Pos)),
			  plotRaw = FALSE,
			  plotWindows = FALSE,
			  method = 4)
           }
	   win <- cbind(Chrom = rep(scaffolds[i], nrow(win$windowData)), win$windowData)
	   write.table(win, file = "genwin.out.csv", append = T, quote = F, sep = ",", row.names = F, col.names = !file.exists("genwin.out.csv"))
	   message(paste0("Finished scaffold ", scaffolds[i]))
   }
}
```
Using the above script, we were able to identify 1,156,611 _de novo_ windows.
_Summary stats of the windows_
```R
# Load in window output
windows <- read.csv("genwin.out.csv", header = T)
   #1,202,291 windows
windows <- subset(windows, MeanY != "NA")
   #1,156,611 windows
windows$lengths <- windows$WindowStop - windows$WindowStart + 1
summary(windows$lengths); sd(windows$lengths)
summary(windows$lengths); sd(windows$lengths)

# Save output
write.table(windows, file = "windows.fst.csv", quote = F, sep = ",", row.names = F, col.names = T)

# Build Plots
library(ggplot2)
library(patchwork)
p1 <- ggplot(windows, aes(x = lengths)) + geom_histogram() + scale_x_continuous(trans='log10')
p2 <- ggplot(windows, aes(x = SNPcount)) + geom_histogram() + scale_x_continuous(trans='log10')
png(file = "windows.hist.pdf", width = 16, height = 10)
p1 + p2
dev.off()
```
<div align="center">
	
| Metric | Min | Mean (SD) | Max |
| --- | --- | --- | --- |
| Window Length | 1 | 977.6 (3057.6) | 299,828 |
| SNP Count | 1 | 12.4 (19.5) | 1,670 |
</div>

<p align="center">
  <img src="figures/windows.hist.png" alt="Histograms" width="750">
</p>
<p align="center">
  <sup>Histograms of the window lengths inferred from _GenWin_ (left) and the number of SNPs per window (right).</sup>
</p>

### Step 2: Find Outlier Windows with Signifcant SNPs
We defined outlier windows based on two criteria:
1. Top 1% of F<sub>ST</sub> values (F<sub>ST</sub> ≥ 0.1733967)
2. Containing at least two significant SNPs
   - as defined with the Fisher's exact test using a false discovery rate (FDR) < 0.01
These are overall very conservative criteria, and consistent with previous studies using poolseq (e.g. [Sly et al. 2022](https://doi.org/10.1073/pnas.2120482119))

_Identify the outlier windows matching the two criteria above_
```R
# Load in window output
windows <- read.csv("windows.fst.csv", header = T)
   #1,156,611 windows

# Get percentiles (top 5%, 2.5%, 1%)
p <- quantile(windows$MeanY, c(0.95, 0.975, 0.99))
#       95%     97.5%       99% 
# 0.1129424 0.1370687 0.1733967 

# Filter for just the top 1% outlier windows
top1.win <- subset(windows, MeanY >= p[3])

# Load in Fisher exact test data
SNPs <- read.table("pools.fet", sep = "\t", header = F)[,c(1,2,6)]
colnames(SNPs) <- c("Chrom", "Pos", "Q")
   #14,439,829 sites

# Calculate FDR and filter for significant SNPs
SNPs$Q <- as.numeric(gsub("1:2=", "", SNPs$Q))
SNPs$P <- 10^(-SNPs$Q) # convert Q value from Popoolation2 to P value
SNPs$FDR <- p.adjust(SNPs$P, method = "fdr") # FDR corrected P-values
SNPs.sig <- subset(SNPs, FDR < 0.01)
   # 1,890 SNPs
```
