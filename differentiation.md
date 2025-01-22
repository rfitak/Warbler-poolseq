# Detecting a genetic effect by season
In this section we tested the hypothesis that there is a larger genetic contribution to the spring migration timing as opposed to the autumn migraton timing.  In essense, we expected the differentiation to be higher between _early_ vs. _late_ chronoptypes in the spring as opposed to the fall. To quantify differentiation, we used both F<sub>ST</sub> and allele frequency differentiation.
The key files needed are:
- _pools.fst_: per site F<sub>ST</sub> calculated with _popoolation2_
- _pools.fet_: per site Fisher's exact test of allele frequency differentiation calculated with _popoolation2_
- _pools_pwc_: per site pairwise allele frequency differences calculated with _popoolation2_

We are only interested in the **1v2** (Early Spring vs Late Spring) and **3v4** (Early Autumn vs Late Autumn) comparisons.

```R
# Read in PoPoolation files
fst <- read.table("pools.fst", header = F, sep = "\t")
fet <- read.table("pools.fet", header = F, sep = "\t")
diff <- read.table("pools_pwc", header = F, sep = "\t")
   #14,439,829 sites

# Fix fst and fet column names
colnames(fst) <- c("Chrom", "Pos", "n", "n2", "n3", "Fst1v2", "Fst1v3", "Fst1v4", "Fst2v3", "Fst2v4", "Fst3v4")
colnames(fet) <- c("Chrom", "Pos", "n", "n2", "n3", "Fet1v2", "Fet1v3", "Fet1v4", "Fet2v3", "Fet2v4", "Fet3v4")
colnames(diff) <- c("Chrom", "Pos", "rc", "allele_count", "allele_states", "delete_sum", "snp_type", "most_variable_allele", "diff1v2", "diff1v3", "diff1v4", "diff2v3", "diff2v4", "diff3v4")

# Remove the "1v2:" from values and make them numeric
fst[,c(6:11)] <- apply(fst[,c(6:11)], 2, function(x) {as.numeric(gsub("[123]:[234]=", "", x))})
fet[,c(6:11)] <- apply(fet[,c(6:11)], 2, function(x) {as.numeric(gsub("[123]:[234]=", "", x))})


# New function for the Standard deviation of a population
sd.pop <- function(x) {sqrt(sum((x - mean(x))^2)/length(x))}

####### Z test for Fst (spring vs Fall)
library(BSDA)
z.test(fst$Fst1v2, fst$Fst3v4, sigma.x = sd.pop(fst$Fst1v2), sigma.y = sd.pop(fst$Fst3v4))

	Two-sample z-Test

data:  fst$Fst1v2 and fst$Fst3v4
z = 118.8, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.002447763 0.002529886
sample estimates:
 mean of x  mean of y 
0.04992974 0.04744091 

z.test(diff$diff1v2, diff$diff3v4, sigma.x = sd.pop(diff$diff1v2), sigma.y = sd.pop(diff$diff3v4))

	Two-sample z-Test

data:  diff$diff1v2 and diff$diff3v4
z = 25.813, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.0009271907 0.0010795622
sample estimates:
mean of x mean of y 
0.1103028 0.1092994 
```

Here is some extra R plotting code:
```R
library(ggplot2)
library(patchwork)
library(reshape2)
data <- melt(cbind(Spring=fst$Fst1v2, Autumn=fst$Fst3v4))
ggplot(data, aes(x = value, fill = as.factor(Var2))) + geom_histogram(position = "identity", alpha = 0.2) + scale_x_continuous(trans='log10')
ggplot() +
   geom_histogram(aes(x = subset(data, Var2 == "Spring", value), fill = "data1"), alpha = 0.2) +
   geom_histogram(aes(x = subset(data, Var2 == "Autumn", value), fill = "data2"), alpha = 0.2) +
   scale_x_continuous(trans='log10') + scale_fill_manual(values = c("data1" = "red", "data2" = "green"))

p1 <- ggplot(data, aes(x = value, y = as.factor(Var2))) + geom_histogram()
data2 <- melt(cbind(Spring=diff$diff1v2, Autumn=diff$diff3v4))
p2 <- ggplot(data2, aes(x = value, y = as.factor(Var2))) + geom_histogram()

# Manhattan Plot
library(qqman)
data <- read.csv("windows.fst.csv", header=T)
data$CHR <- as.numeric(factor(data$Chrom, labels = 1:length(unique(data$Chrom))))
data$BP <- (data$WindowStop + data$WindowStart) / 2
data$SNP = paste0(data$Chrom, "-", data$BP)
data.sig <- read.csv("outlier-windows.fst.csv", header = T)
data.sig$BP <- (data.sig$WindowStop + data.sig$WindowStart) / 2
data.sig$SNP = paste0(data.sig$Chrom, "-", data.sig$BP)
pdf(file = "Rplots2.pdf", width = 16, height = 10)
manhattan(data,
	chrlabs = rep("", length(unique(data$CHR))),
	p = "MeanY",
	logp = F,
	cex = 1.0,
	col = c("blue4", "orange3"),
	suggestiveline = F,
	genomewideline = F,
	xlab = "Scaffold",
	ylab = "FST",
	xaxt = "n",
	highlight = data.sig$SNP)
dev.off()
```
