# PrimeCUTR: Neopeptide Prediction from Somatic Mutations
PrimeCUTR is an R package designed for predicting neopeptides from somatic mutation calls in cancer genomic sequencing data (Sng et al., 2024). The package includes functionality to handle start-gain, stop-loss, missense, and frameshift mutations. It requires Variant Effect Predictor (VEP) annotated variant call format (VCF) files as input and produces neopeptide sequences suitable for downstream analysis, such as MHC binding prediction using tools like netMHC.

## Installation
To install PrimeCUTR, you can use the following commands in R:

```R
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install PrimeCUTR from GitHub
devtools::install_github("christophersng/primeCUTR")
```

## Getting started
To use PrimeCUTR, load the package in your R script or R Markdown document:

```R
library(primeCUTR)
```
<br>

The key function in PrimeCUTR is **_get.peptide_** which can be used interactively in R to obtain the neopeptide and metadata for a given transcript and mutation combination. 

```R
get.peptide("ENST00000539214","c.-61C>T",build = 38,check_startgains = TRUE)
```
<br>

In order to process an entire .vcf file of somatic mutations, use **_get.orfs_**:

```R
#two example VCF files are bundled with this package as can be accessed like so:
input_vcf_path <- primeCUTR_example("vep_HCC1937.vcf")
#OR
input_vcf_path <- primeCUTR_example("vep_MCF-7.vcf")

#Run the function
get.orfs(input_vcf_path,"./output_dir/",build=38)
```
<br>
A wrapper R script can be found here (https://github.com/christophersng/primeCUTR_scripts) for simple loop batch running of multiple .vcf files in a folder.

Note that the VCF files for input to PrimeCUTR are dependent on VEP-annotation 
(https://www.ensembl.org/info/docs/tools/vep/index.html) with the `--hgvs` flag option selected. For example:

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

For small batches of somatic mutation data, Ensembl VEP online is freely accessible to all users at: https://www.ensembl.org/Tools/VEP.
Simply ensure that the `HGVS` option is checked.

## Outputs
**_get.peptide_** used interactively returns the following outputs in the R console:<br>
(Run `?get.peptide` for more information on input arguments and usage of the function.)<br>
<br>
Core

 * `$peptide.seq`: The full mutated neopeptide sequence.
 * `$start`: The index position for the first mutated/frameshifted peptide.
 * `$end`: The index position for the last mutated/frameshifted peptide.
 * `$peptides_msg`: Any additional metadata on the nature of the neoORF, and/or error messages.

Optional (for start-gains only), if `check_startgain=T`

 * `$wt_kz`: Native coding start site Kozak sequence score<sup>✝</sup>.
 * `$mut_kz`: Mutant start-gain Kozak sequence score<sup>✝</sup>.
 * `$relative_strength_overlap`: The difference in Kozak sequence score between the neoORF and any overlapping uORF.
 * `$atg_count`: The number of ATG uORFs in the wild-type 5'UTR.
 * `$overlap_present`: If `TRUE` indicates the presence of an overlapping wild-type uORF in the same reading frame.
 * `$overlap_pep`: The overlapping sequence.
 * `$perc_overlap`: How much of the neoORF consists of overlap.
 * `$relative_strength_native`: The difference in Kozak sequence score between the neoORF and the native coding start site.
 * `$cds_overlap_expected`:  If `TRUE` suggests in-frame overlap with the coding sequence of another isoform.
 * `$overlap_cds_transcripts`: Comma-delimited list of overlapping isoforms.

<br>
<sup>✝</sup>Calculated per the schema in Whiffin, Nicola, et al. "Characterising the loss-of-function impact of 5’untranslated region variants in 15,708 individuals." Nature communications 11.1 (2020): 2523.
<br>
<br>

**_get.orfs_** produces 3 output folders in the designated `output` folder.<br>
(Run `?get.peptide` for more information on input arguments and usage of the function.)

```
output
  ├──log
  │   └──study_id_patient_id_YYYY_mmm_dd_HH_mm_ss.txt
  ├──netMHC_inputs
  │   ├──study_id_patient_id_missenseclass_9mer_peptides.YYYY_mmm_dd_HH_mm_ss.txt
  │   ├──study_id_patient_id_missenseclass_10mer_peptides.YYYY_mmm_dd_HH_mm_ss.txt
  │   ├──study_id_patient_id_missenseclass_11mer_peptides.YYYY_mmm_dd_HH_mm_ss.txt
  │   ├──study_id_patient_id_frameshiftclass_9mer_peptides.YYYY_mmm_dd_HH_mm_ss.txt
  │   └── ...
  └──orfs
      ├──study_id_patient_id_missenseORFs_YYYY_mmm_dd_HH_mm_ss.tsv
      ├──study_id_patient_id_frameshiftORFs_YYYY_mmm_dd_HH_mm_ss.tsv
      ├──study_id_patient_id_startgainORFs_YYYY_mmm_dd_HH_mm_ss.tsv
      └──study_id_patient_id_stoplossORFs_YYYY_mmm_dd_HH_mm_ss.tsv
```

The `netMHC_input` files are FASTA-format text files containing 9, 10 and 11-mers which include at least one mutated/frameshifted residue. 
These can be used for downstream MHC binding prediction. 
<br><br>
The `orfs` files are tab-separated file containing neo-peptides per mutation class per sample. Apart from those outputs in `get.peptide`,
the tab-separated files also contain the following additional columns:

 * `$key`: patient_id--hgnc_symbol--ensembl_transcript_id--chr:pos::ref:alt--hgvs_simple_annotation--mutant-sequence
 * `$hgnc_symbol`: HUGO Gene Nomenclature Committee Symbol
 * `$ensembl_transcript_id`: Ensembl Transcript ID as provided by the VEP annotation
 * `$peptide`: The sequence of mutated/frameshifted amino acid residues
 * `$peptide.tryptic.ends`: The sequence of mutated/frameshifted amino acid residues with additional normal flanking residues, up to the next trypsin cleavage site with one missed cleavage per side
 * `$peptide.normal.ends`: The sequence of mutated/frameshifted amino acid residues with 10 additional normal flanking residues

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
Contributions and feedback on PrimeCUTR are welcome. If you encounter issues or have suggestions for improvements, please open an issue.

## Citation
Sng CCT, Kallor AA, Simpson BS, Bedran G, Alfaro J, Litchfield K. Untranslated regions (UTRs) are a potential novel source of neoantigens for personalised immunotherapy. _Frontiers in Immunology_. 2024 Mar 15;15:1347542.
https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1347542/full
