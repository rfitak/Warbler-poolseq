# Alternative poolseq pipeline with _Grenedalf_

After completing the traditional *popoolation2* pipeline above, A quick literature search revealed that there is a new software for poolseq analysis, [grenedalf](https://github.com/lczech/grenedalf) ([Czech et al. 2024](https://doi.org/10.1093/bioinformatics/btae508)). This new software was developed specifically in conjunction with Robert Kofler, the author of both *popoolation* and *popoolation2*. The paper is a great read, as it is a thorough treatment of the population genetic parametes used in the sofware (i.e., $\pi$, F<sub>ST</sub>, Tajima's D, etc) and how their derivation needs to be modified for poolseq data. <u>It turns out that the Tajima's D cannot actually be properly derived for poolseq data, and any study published previously using this metric from pooled data must consider a re-analysis</u>:
> "In addition, we noticed several implementation bugs in POPOOLATION (Kofler *et al.* 2011a), up until and including v1.2.2 of the tool. We discussed these with the authors, and the bugs have since been fixed (pers. comm. with R. Kofler). If conclusions of studies depend on numerical values of Tajima’s *D*computed with POPOOLATION, we recommend reanalyzing the data. Due to these statistical issues, we generally advise to be cautious when applying and interpreting Tajima’s *D* with Pool-seq data."

Although our purple martin analysis is not using metrics of variation, and the F<sub>ST</sub> implementation is generally unchanged, I still wanted to run this new pipeline on our pools as well. *grenedalf v0.6.2* was pretty easy to install and run, and actually much more user friendly and scalable than *popoolation2*. It also has an improved interface for defining windows and their denomination to decrease biases (read their paper and supplementary methods). I highly recommend it. This code is below.

```bash
# Max coverage = 2.5X the mean depth: (19.2 + 22.3 + 25.5 + 22.6) / 4 * 2.5 = 56

# Get allele frequencies
grenedalf \
	frequency \
	--sam-path B1.clean.sorted.bam \
	--sam-path B2.clean.sorted.bam \
	--sam-path C1.clean.sorted.bam \
	--sam-path C2.clean.sorted.bam \
	--sam-min-map-qual 20 \
	--sam-min-base-qual 20 \
	--reference-genome-fasta Aarun.fa \
	--write-sample-counts \
	--write-sample-read-depth \
	--write-sample-ref-freq \
	--separator-char comma \
	--file-prefix grenedalf. \
	--verbose \
	--threads 32

# Get Fst metrics
grenedalf \
        fst \
	--sam-path B1.clean.sorted.bam \
        --sam-path B2.clean.sorted.bam \
        --sam-path C1.clean.sorted.bam \
        --sam-path C2.clean.sorted.bam \
	--sam-min-map-qual 20 \
        --sam-min-base-qual 20 \
        --reference-genome-fasta Aarun.fa \
        --filter-sample-min-count 2 \
        --filter-sample-max-count 0 \
        --filter-sample-min-read-depth 5 \
        --filter-sample-max-read-depth 56 \
        --filter-total-only-biallelic-snps \
        --window-type interval \
        --window-interval-width 25000 \
        --window-interval-stride 25000 \
        --window-average-policy valid-loci \
        --method unbiased-nei \
        --pool-sizes 18 \
        --file-prefix grenedalf. \
        --verbose \
        --threads 32
```

_Parameters Explained:_

- ***--sam-path file*** :: path to BAM file of mapped reads, cleaned of duplicates, poorly mapepd reads, etc.
- ***--sam-min-map-qual 20*** :: minimum mapping quality score to consider a read
- ***--sam-min-base-qual 20*** :: minimum base quality score to consider
- ***--reference-genome-fasta Aarun.fa*** :: fasta file of the reference genome
- ***--write-sample-counts*** :: write 'REF_CNT' and 'ALT_CNT' columns per sample, containing the REF and ALT base counts at the position for each sample
- ***--write-sample-read-depth*** :: write a 'DEPTH' column per sample, containing the read depth (sum of REF and ALT) counts of each sample.
- ***--write-sample-ref-freq*** :: write a 'FREQ' column per sample, containing the reference allele frequency, computed as REF/(REF+ALT) of the counts of each sample
- ***--separator-char comma*** :: Separator char between fields of output tabular data {comma, tab, space, semicolon}
- ***--filter-sample-min-count 2*** :: Minimum base count for a nucleotide (in `ACGT`) to be considered as an allele
- ***--filter-sample-max-count 0*** :: Maximum base count for a nucleotide (in `ACGT`) to be considered as an allele
- ***--filter-sample-min-read-depth 5*** :: Minimum read depth expected for a position in a sample to be considered covered. If the sum of nucleotide counts (in `ACGT`) at a given position in a sample is less than the provided value, the sample is ignored at this position.
- ***--filter-sample-max-read-depth 56*** :: Maximum read depth expected for a position in a sample to be considered covered. If the sum of nucleotide counts (in `ACGT`) at a given position in a sample is greater than the provided value, the sample is ignored at this position.
- ***--filter-total-only-biallelic-snps*** :: Filter out any positions that do not have exactly two alleles across all samples. That is, after applying all previous filters, if not exactly two counts (in `ACGT`) are non-zero in total across all samples, the position is not considered a biallelic SNP, and ignored.
- ***--window-type interval*** :: Type of window to use. Depending on the type, additional options might need to be provided.
  - `interval`: Typical sliding window over intervals of fixed length (in bases) along the genome
  - `queue`: Sliding window, but instead of using a fixed length of bases along the genome, it uses a fixed number of positions of the input data. Typically used for windowing over variant positions such as (biallelic) SNPs, and useful for example when SNPs are sparse in the genome.
  - `single`: Treat each position of the input data as an individual window of size 1. This is typically used to process single SNPs, and equivalent to `interval` or `queue` with a width/count of 1, except that positions that are removed by some filter are skipped.
  - `regions`: Windows corresponding to some regions of interest, such as genes. The regions need to be provided, and can be overlapping or nested as needed.
  - `chromosomes`: Each window covers a whole chromosome.
  - `genome`: The whole genome is treated as a single window.
- ***--window-interval-width 25000*** :: Required when using `--window-type interval`: Width of each window along the chromosome, in bases
- ***--window-interval-stride 25000*** :: When using `--window-type interval`: Stride between windows along the chromosome, that is how far to move to get to the next window
- ***--window-average-policy valid-loci*** :: Denominator to use when computing the average of a metric in a window:
  - `window-length`: Simply use the window length, which likely underestimates the metric, in particular in regions with low coverage and high missing data.
  - `available-loci`: Use the number of positions for which there was data at all (that is, absent or missing data is excluded), independent of all other filter settings.
  - `valid-loci`: Use the number of positions that passed all quality and numerical filters (that is, excluding the SNP-related filters). This uses all positions of high quality, and is the recommended policy when the input contains data for all positions.
  - `valid-snps`: Use the number of SNPs only. This might overestimate the metric, but can be useful when the data only consists of SNPs.
  - `sum`: Simply report the sum of the per-site values, with no averaging applied to it. This can be used to apply custom averaging later.
  - `provided-loci`: Use the exact loci provided via `--window-average-loci-bed` or `--window-average-loci-fasta` to determine the denominator for the window averaging, by counting all positions set in this mask in the given window.
- ***--method unbiased-nei*** :: FST method to use for the computation:
  - `unbiased-nei` or `unbiased-hudson` The unbiased pool-sequencing statistic (in two variants, following the definition of Nei, and the definition of Hudson)
  - `kofler` the statistic by Kofler et al of PoPoolation2
  - `karlsson` the asymptotically unbiased estimator of Karlsson et al (which is also implemented in *PoPoolation2*)
  - All except for the Karlsson method also require `--pool-sizes` to be provided.
- ***--pool-sizes pool.sizes*** :: Pool sizes for all samples that are used (not filtered out). These are the number of haploids, so 100 diploid individuals correspond to a pool size of 200. Either:
  - a single pool size that is used for all samples, specified on the command line, or
  - a path to a file that contains a comma- or tab-separated list of sample names and pool sizes, with one name/size pair per line, in any order of lines.
- **--file-prefix grenedalf.**. :: File prefix for output files
- **--verbose** :: use verbose output reporting
- **--threads 32** :: number of parallel threads to use
