# Pool-seq Analysis with _Popoolation2_
Now that all our pooled BAM files are prepared, we are ready to run the pool-seq analysis pipeline.  We used [popoolation2 v1.201](https://sourceforge.net/p/popoolation2/wiki/Main/) ([Kofler et al. 2011](https://doi.org/10.1093/bioinformatics/btr589)). For additional details, please see their review paper in *Nature Reviews Genetics* ([Schlötterer et al. 2014](https://doi.org/10.1093/bioinformatics/btr589)) or this useful tutorial by Joanna Griffiths in the ***MarineOmics*** group (https://marineomics.github.io/POP_03_poolseq.html).

Below, the mandatory parameter `max-coverage` was calculated as 2.5X the coverage per pool.
-  19.2 | 22.3 | 25.5 | 22.6 x 2.5 = 48 | 56 | 64 | 57

```bash
# Build the pileup of the 4 pools
samtools mpileup -B \
   B1.clean.sorted.bam \
   B2.clean.sorted.bam \
   C1.clean.sorted.bam \
   C2.clean.sorted.bam \
   > pools.mpileup

# Generate synchronized file
   # this is the multithreaded, ultrafast version only found in v1.201 of popoolation2
java -ea -Xmx7g -jar popoolation2_1201/mpileup2sync.jar \
   --input pools.mpileup \
   --output pools.sync \
   --fastq-type sanger \
   --min-qual 20 \
   --threads 32

# Calculate allele frequencies and differences
perl popoolation2_1201/snp-frequency-diff.pl \
   --input pools.sync \
   --output-prefix pools \
   --min-count 2 \
   --min-coverage 5 \
   --max-coverage 48:56:64:57

# Calculate pairwise Fst values
perl popoolation2_1201/fst-sliding.pl \
   --input pools.sync \
   --output pools.fst \
   --min-count 2 \
   --min-coverage 5 \
   --max-coverage 48:56:64:57 \
   --min-covered-fraction 1 \
   --window-size 1 \
   --step-size 1 \
   --pool-size 18 \
   --suppress-noninformative

# Calculate exact test of allele frequency differences
perl popoolation2_1201/fisher-test.pl \
   --input pools.sync \
   --output pools.fet \
   --min-count 2 \
   --min-coverage 5 \
   --max-coverage 48:56:64:57 \
   --min-covered-fraction 1 \
   --window-size 1 \
   --step-size 1 \
   --window-summary-method multiply \
   --suppress-noninformative
```

_Parameters Explained:_
- ***--fastq-type sanger*** :: quality scores are encoded in the standard Sanger format
- ***--min-qual 20*** :: minimum base quality score to consider
- ***--threads 32*** :: use 32 threads (multithreaded)
- ***--min-count 2*** :: the minimum count of the minor allele. used for SNP identification. SNPs will be identified considering all populations simultanously.
- ***--min-coverage 5*** :: the minimum coverage; used for SNP identification, the coverage in ALL populations has to be higher or equal to this threshold, otherwise no SNP will be called.
- ***--max-coverage 57*** :: The maximum coverage; All populations are required to have coverages lower or equal than the maximum coverage; Mandatory; The maximum coverage may be provided as one of the following:
  -  '500' :: a maximum coverage of 500 will be used for all populations
  -  '300,400,500' :: a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on...
  -  '2%' :: the 2% highest coverages will be ignored, this value is independently estimated for every population
- ***--min-covered-fraction 1*** :: the minimum fraction of a window being between min-coverage and max-coverage in ALL populations
- ***--window-size 1*** :: the size of the sliding window. Measured in "--window-unit"
- ***--step-size 1*** :: the size of the sliding window steps. Measured in "--window-unit"
- ***--pool-size 18*** :: the size of the population pools; May be provided for each population individually; mandatory parameter. NOTE THAT THIS IS THE # OF CHROMOSOMES (2 x # individuals for diploids)
  - --pool-size 500 .. all populations have a pool size of 500
  - --pool-size 500:300:600 .. first population has a pool size of 500, the seccond of 300 etc; the number of pools has to fit with the number of populations provided in the file
- ***--window-summary-method multiply*** :: Specify the method by which the p-values of individual SNPs should be summarized; possible: geometricmean | median | multiply
- ***--suppress-noninformative*** :: Suppress output for windows with no SNPs or insufficient coverage
