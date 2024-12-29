<h3><p align="center">Supplementary Methods for:</p></h3>
<h2><p align="center">The influence of genes on the migratory timing of great reed warblers (<i>Acrocephalus arundinaceus</i>)</p></h2>

<I><h5>Emily R. Fackler<sup>1</sup>, Dmitry Kishkinev<sup>2</sup>, Petr Prochazka<sup>3</sup>, Robert R. Fitak<sup>2</sup></h5></I>

1. Department of Biology, Genomics and Bioinformatics Cluster, University of Central Florida, Orlando, FL 32816
2. Keele University, United Kingdom
3. Institute of Vertebrate Biology, Czech Republic

<br>
<p align="center">
  <img src="images/xxxxxx.jpg" alt="title" width="500">
</p>
<p align="center"><sup>Figure caption.</sup>
</p>
<br>

## Raw Data
All project, sample information, and raw sequencing data are publicly available and can be accessed through the NCBI BioProject database (accession [PRJNA995180](https://www.ncbi.nlm.nih.gov/bioproject/995180)).

| Pool | Abbreviation | BioSample Accession | SRA Accession |
| --- | ---| --- | --- |
| Early Spring | B1 | [SAMN37026500](https://www.ncbi.nlm.nih.gov/biosample/37026500) | [SRR25778514](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR25778514) |
| Late Spring | B2 | [SAMN37026501](https://www.ncbi.nlm.nih.gov/biosample/37026501) | [SRR25778513](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR25778513) |
| Early Autumn | C1 | [SAMN40709721](https://www.ncbi.nlm.nih.gov/biosample/40709721) | [SRR28726269](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR28726269) |
| Late Autumn | C2 | [SAMN40709722](https://www.ncbi.nlm.nih.gov/biosample/40709722) | [SRR28726268](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR28726268) |


***
___This GitHub repository contains a summary of the various code, software, and data analysis pipelines used for the aforementioned study of the great reed warbler. The contents represented here are only to be used as an example and not intended to be comeprehensive. The authors make no representation about the suitability or accuracy of this code, software, or data for any purpose, and make no warranties, either expressed or implied, for a particular purpose or that the use of this software or data will not infringe any third party patents, copyrights, trademarks, or other rights. The code, software and data are provided "as is". All content is hereby registered under the GNU General Public License v3.0, see [LICENSE](./LICENSE). Any publication that significantly relies upon the use of the content generated herein shall appropriately cite:___

<p align="center">Fackler et al. (in prep.) The influence of genes on the migratory timing of great reed warblers (<i>Acrocephalus arundinaceus</i>). TBD</p>

***
  
<h2><p align="center">Table of Contents</p></h2>
<div align="center">
 
[Read processing and cleaning](./read_processing.md)

[Pool-seq analysis: <i>Popoolation2</i>](./popoolation2.md)

[Pool-seq analysis: <i>grenedalf</i>](./grenedalf.md)

[Window analysis: <i>GenWin</i>](./windows.md)

[Enrichment analysis](./enrichment.md)

[Genetic Effect by Season](./differentiation.md)

[Alternative Reference Genome Complete Analysis](./alternate-genome.md)

</div>

***

<h2><p align="center">Summary of the Programs/Software Used in this Study</p></h2>  

| Program | Version | Citation |
| --- | --- | --- |
| FASTP | 0.20.0 | [Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890.](https://doi.org/10.1093/bioinformatics/bty560) |
| SAMTOOLS | 1.18 | [Li, H. et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079.](https://doi.org/10.1093/bioinformatics/btp352) |
| BWA | 0.7.17 | [Li and Durbin (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754–1760.](https://doi.org/10.1093/bioinformatics/btp324) |
| POPOOLATION2 | 1.201 | [Kofler et al. (2011). PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27, 3435–3436.](https://doi.org/10.1093/bioinformatics/btr589) |
| GRENEDALF | 0.6.2 | [Czech et al. (2024). grenedalf: population genetic statistics for the next generation of pool sequencing. Bioinformatics 40, btae508.](https://doi.org/10.1093/bioinformatics/btae508) |
| GENWIN | 1.0 | [Beissinger et al. (2015). Defining window-boundaries for genomic analyses using smoothing spline techniques. Genetics Selection Evolution 47, 30.](https://doi.org/10.1186/s12711-015-0105-9) |
| BEDTOOLS | 2.27.1 | [Quinlan and Hall (2010)., BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841–842.](https://doi.org/10.1093/bioinformatics/btq033) |

