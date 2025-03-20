#!/bin/bash

#SBATCH --job-name=plink.clump
#SBATCH --nodes=1
#SBATCH --ntasks=1
#sbatcH --cpus-per-task=5
#sbATCH --mem=30G
#SBATCH --time=5:00:00
#SBATCH --partition=general
#SBATCH --output .../gwas_threshold/log_out/%x.%j.out
#SBATCH --error .../gwas_threshold/log_out/%x.%j.err
#SBATCH --array=1-25 # number of metabolites to run in the list

## code USAGE:
#for chr in {1..22}; do sbatch --export=chr=$chr .../scripts/03a_per_trait_metabolites_plink_clump.sbatch; done

## Process input data by task
task_id=$SLURM_ARRAY_TASK_ID
chr=$chr

echo "plink clump on trait " $task_id "and Chr " $chr

## set-up parameters for plink
gwas_cutoff='8.87e-9'
maf_cutoff=1e-3
analysis=additive

## Contain list of trait names to run plink & summary stats
# metabolites  traits
trait_lst=".../02_selected_metabolites_irnt_bothsex.in"
gwas_stats=".../neale_metabolites_gwas_stats"

## other traits 
#trait_lst=".../01_selected_gwas_48_irnt_bothsex_v3.in"
#gwas_stats=.../UKB_neale_R2/results_uncompress

## other directories
base_dir="..."
ukb_geno_dir="..."
ukb_10k_ld="..."
gwas_files="..."

## output dir
output_path=$base_dir/out_plinkclump_maf${maf_cutoff}_${analysis}_gwasthresh${gwas_cutoff}

## create work dir
workdir=`mktemp -d $base_dir/work_dir/tmpXXXXXXX`
echo "Temporary directory to create is" $workdir
mkdir -p $output_path

## load GCTA
module load gcta

## extract the traitName from the filename (remove the .tsv extension)
filename=$(awk -v trait=$task_id '$1==trait{print $2}' $trait_lst )
traitName=$(basename "$filename" .tsv)

## cp files to the work dir
cp "$gwas_stats/$filename" $workdir
## Go to work dir
cd $workdir

## modify SNP names and prune based on MAF
awk 'NR==FNR {a[$2]=$0; next} $1 in a {print a[$1], $0}' \
    <(cat "$filename" | awk '{print $2, $3, $8, $9, $11, $5, $1}' | \
            sed 's/:/ /g' | awk '{print $7, $7":"$8, $8, $1, $9, $6, $3, $4, $5}' | sed '1d') \
    <(awk '{print $1 ":" $4, $0}' "$ukb_10k_ld/ukbEURu_imp_chr${chr}_v3_impQC_10k_mac1.bim") | \
            awk '{print $1, $12, $3, $4, $5, $6, $7, $8, $9}' | \
    sed '1i CHR SNP BP A1 A2 NMISS BETA SE P' > "${traitName}.ma"

## filter GWAS stats based on MAF
awk '$4>1e-3{print $0}' ${traitName}.ma > ${traitName}.ma.sub

echo "plink clump on trait " $filename

## run plink clumping
plink --bfile $ukb_10k_ld/ukbEURu_imp_chr${chr}_v3_impQC_10k_mac1 \
    --clump ${traitName}.ma.sub \
    --clump-p1 ${gwas_cutoff} \
    --clump-r2 0.01 \
    --clump-kb 250  \
    --out ${traitName}.chr${chr}

## cp results and delete dir.
rsync -av *.clumped $output_path  && echo "Files copied to: $output_path"

## delete work_dir
cd ..
[ -d "$workdir" ] && rm -rf "$workdir"
