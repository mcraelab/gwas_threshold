#!/bin/bash

#SBATCH --job-name=cojo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#sbatcH --cpus-per-task=5
#sbATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --partition=general
#SBATCH --output .../log_out/%x.%j.out
#SBATCH --error .../log_out/%x.%j.err
#SBATCH --array=1-25

## USAGE:
#for chr in {1..22}; do sbatch --export=chr=$chr .../scripts/03a_cojo3_per_trait_re_analysis.sbatch; done

## Process input data by task
task_id=$SLURM_ARRAY_TASK_ID
chr=$chr

echo "COJO on trait " $task_id "and Chr " $chr

## set-up parameters for cojo
gwas_cutoff='8.87e-9'
maf_cutoff=1e-3
analysis=additive

## Contain list of trait names to run cojo & summary stats
# metabolites
trait_lst=".../02_selected_metabolites_irnt_bothsex.in"
gwas_stats="..."

## other traits 
#trait_lst=".../01_selected_gwas_48_irnt_bothsex_v3.in"
#gwas_stats="..."

## Other directories
base_dir="..."
ukb_geno_dir="..."
ukb_10k_ld="..."
gwas_files="..."

## output dir
#output_path=$base_dir/out_cojo_maf${maf_cutoff}_gwasthresh${gwas_cutoff} # for 5.0e-8 cut-off
output_path=$base_dir/out_cojo_maf${maf_cutoff}_${analysis}_gwasthresh${gwas_cutoff}
## create work dir
workdir=`mktemp -d $base_dir/work_dir/tmpXXXXXX`
echo "Temporary directory to create is" $workdir
mkdir -p $output_path

## load GCTA
module load gcta
## extract the traitName from the filename (remove the .tsv extension)
filename=$(awk -v trait=$task_id '$1==trait{print $2}' $trait_lst )
traitName=$(basename "$filename" .tsv)

## check and exit
#[[ -f $output_path/${traitName}.chr${chr}.jma.cojo ]] && exit
if [[ -f $output_path/${traitName}.chr${chr}.jma.cojo ]]; then
    echo "The file: " ${traitName}.chr${chr}.jma.cojo " exists...exiting"
    exit
fi

## cp files to the work dir
cp "$gwas_stats/$filename" $workdir
## Go to work dir
cd $workdir

## prune data
awk 'NR==FNR {a[$1]=$0; next} $1 in a {print a[$1], $0}' <(cat "$filename" | \
    awk '{print $2, $3, $8, $9, $11, $5, $1}' | sed 's/:/ /g' | \
    awk '{print $7 ":" $8, $1, $9, $2, $3, $4, $5, $6}' | \
    sed '1d' | \
    sed '1i SNP minor_allele A2 minor_AF beta se pval n_complete_samples')  \
     <(awk '{print $1 ":" $4, $0}' "$ukb_10k_ld/ukbEURu_imp_chr${chr}_v3_impQC_10k_mac1.bim") | \
     awk '{print $11, $2, $3, $4, $5, $6, $7, $8}' > "${traitName}.ma"
## filter GWAS stats based on MAF
awk '$4>1e-3{print $0}' ${traitName}.ma > ${traitName}.ma.sub

echo "Cojo on trait " $filename
## run COJO
gcta64 --bfile $ukb_10k_ld/ukbEURu_imp_chr${chr}_v3_impQC_10k_mac1 \
       --maf 0.001 \
       --cojo-file ${traitName}.ma.sub \
       --cojo-slct \
       --thread-num 5 \
       --cojo-wind 1000 \
       --cojo-p $gwas_cutoff \
       --out ${traitName}.chr${chr}
## cp results and delete dir
rsync -av *.jma.cojo $output_path
cd ..
rm -rf $workdir
