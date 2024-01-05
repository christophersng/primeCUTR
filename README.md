# PrimeCUTR: Neopeptide Prediction from Somatic Mutations
PrimeCUTR is an R package designed for predicting neopeptides from somatic mutation calls in cancer genomic sequencing data (Sng et al., 2024 - under peer-review). The package includes functionality for identifying start-gain, stop-loss, missense, and frameshift mutations. It requires Variant Effect Predictor (VEP) annotated variant call format (VCF) files as input and produces neopeptide sequences suitable for downstream analysis, such as MHC binding prediction using tools like netMHC.

## Installation
To install PrimeCUTR, you can use the following commands in R:

```R
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install PrimeCUTR from GitHub
devtools::install_github("christophersng/primeCUTR")
```

## Usage
To use PrimeCUTR, load the package in your R script or R Markdown document:

```R
library(PrimeCUTR)
```

The key function in PrimeCUTR is [_get.peptide_](./man/get.peptide.Rd) which can be used interactively in R to obtain the neopeptide and metadata for a given transcript and mutation combination. 

```R
get.peptide("ENST00000539214","c.-61C>T",build = 38,check_startgains = TRUE)
```

In order to process an entire .vcf file of somatic mutations, use [_get.orfs_](./man/get.orfs.Rd):

```R
#two example VCF files are bundled with this package as can be accessed like so:
input_vcf <- primeCUTR_example("vep_HCC1937.vcf")
get.orfs(input_vcf,"./output_dir/",build=38)

input_vcf <- primeCUTR_example("vep_MCF-7.vcf")
get.orfs(input_vcf,"./output_dir/",build=38)
```

Note that the VCF files for input to PrimeCUTR are dependent on VEP-annotation with the `--hgvs` flag option selected. For example:

```Bash
vep \
--vcf --offline --format vcf \
--af --appris --biotype --buffer_size 5000 --check_existing --distance 5000 \
--hgvs --mane --polyphen b \
--pubmed --regulatory --sift b --species homo_sapiens \
--symbol --transcript_version --tsl --cache \
--dir_cache $cache_loc \
--cache_version $version \
--assembly $assembly \
--fasta $fasta_loc \
--input_file $vcf \
--output_file $output_file
```

## Dependencies
PrimeCUTR relies on the following R packages:

 * BSgenome.Hsapiens.UCSC.hg19
 * BSgenome.Hsapiens.UCSC.hg38
 * dplyr
 * stringr
 * Biostrings

These packages are automatically installed during the PrimeCUTR installation process.
For manual installation of BSgenome packages please follow installation instructions
on the Bioconductor page: https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html

## License
PrimeCUTR is licensed under GPL (>= 3). See the LICENSE file for more details.

## Contribution
Contributions to PrimeCUTR are welcome! If you encounter issues or have suggestions for improvements, please open an issue or submit a pull request.
