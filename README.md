# GWAS significance thresholds simulation

Code for:
GWAS significance thresholds in large cohorts
Evans K Cheruiyot, Tingyan Yang, Allan F McRae
https://doi.org/10.1101/2024.12.09.627629

## Scripts

### analysis/

Bash scripts and slurm jobs for running additive and dominance GWAS, extracting minimum p-values from output, and running GCTA COJO and PLINK clumping

### visualisation

R scripts for plotting figures used in the manuscript

## Notes on usage

 * Directory paths have been replaced with `"..."`.  These need replaced with appropriate values for your system

 * GWAS were run on simluated or permuated phenotypes using a European subset of the UK Biobank cohort. Details on access requirements can be found at https://www.ukbiobank.ac.uk/

 * GWAS results for the GCTA COJO and PLINK clumping analysis were obtained from the Neale Lab (https://www.nealelab.is/uk-biobank).

