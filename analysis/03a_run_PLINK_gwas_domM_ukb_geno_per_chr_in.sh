#!/bin/bash
base_dir="..."
ukb_geno_dir="..."
output_path="..."

[[ -f $output_path/o_chr@chr@_@file@.PHENO1.glm.linear ]] && exit

# Run GWAS analysis for the current chromosome
plink2 --bfile $ukb_geno_dir/ukbEURu_imp_chr@chr@_v3_impQC \
       --pheno $base_dir/pheno/p@file@.pheno \
       --glm dominant allow-no-covars \
       --out o_chr@chr@_@file@ \
       --threads 5
