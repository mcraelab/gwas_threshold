#!/bin/bash --login
#SBATCH --job-name=GWAS-threshold
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --array=1-1000

cd $workdir
file=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
./gcta --mbfile geno_chrs.txt --fastGWA-lr --thread-num 10 --pheno pheno_bmi/p${file}.pheno --out out_additive_bmi/o${file}.out

