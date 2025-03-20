#!/bin/bash

#SBATCH --job-name=domMthresh-GWAS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#sbATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH --output .../UKBGWAS/log_out/%x.%j.out
#SBATCH --error .../UKBGWAS/log_out/%x.%j.err
#SBATCH --array=1-1000

# How to run this batch job:
#for chr in {1..22}; do sbatch --export=chr=$chr .../UKBGWAS/gwas_script/003a_run_gwas_domM_ukb_geno_per_chr.sbatch; done

file=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
echo $file

# run this line to initiate interactive session
#salloc --nodes=1 --cpus-per-task=1 --mem=100G --job-name=TEST --time=10:00:00 --partition=general --account=a_mcrae

# set-up directories
ukb_geno_dir="..."
base_dir="..."
work_dir="..."
output_path="..."
workdir=`mktemp -d $work_dir/tmpXXXXXX`
echo "Temporary directory to create is" $workdir
mkdir -p $output_path

cd $workdir

echo $chr

# Run GWAS analysis for the current chromosome
sed "s/@chr@/$chr/g; s/@file@/$file/g" $base_dir/gwas_scripts/003a_run_gwas_domM_ukb_geno_per_chr_in.sh > run_in.sh

# module load plink2
module load plink

# Run GWAS analysis for the current chromosome & pheno
bash run_in.sh

# combine all chromosomes gwas
#awk 'FNR==1 && NR!=1 { while(/^#CHROM/) getline}1 {print}' *glm.linear > p${file}.glm.linear

rsync -av *.glm.linear $output_path

cd ..
rm -rf $workdir




