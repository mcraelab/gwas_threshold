#!/bin/bash
  
#SBATCH --job-name=AddMinP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#sbATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --output .../log.out/%x.%j.out
#SBATCH --error .../log.out/%x.%j.err
#SBATCH --array=1-1000

file=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
#file=0433

work_dir="..."
output_path=".../out_gwas_min_p_value"
maf_path=".../out_MAF_ukb"

# set-up variables
input_path="..."
analysis=additive_bmi

echo $file
workdir=`mktemp -d $work_dir/tmpXXXXXX`
echo "Temporary directory to create is" $workdir
mkdir -p $output_path

cd $workdir

# Step 1: Filter GWA output based on MAF for each simulated pheno
# MAF 0.1%
awk  'NR==FNR{a[$2]++;next} !($2 in a)' ${maf_path}/MAF_01percent.txt ${input_path}/o${file}.out.fastGWA > o${file}_${analysis}_filter_01pMAF.txt
# MAF 0.5%
awk  'NR==FNR{a[$2]++;next} !($2 in a)' ${maf_path}/MAF_05percent.txt ${input_path}/o${file}.out.fastGWA > o${file}_${analysis}_filter_05pMAF.txt
# MAF 1%
awk  'NR==FNR{a[$2]++;next} !($2 in a)' ${maf_path}/MAF_1percent.txt ${input_path}/o${file}.out.fastGWA > o${file}_${analysis}_filter_1pMAF.txt
# MAF 5% 
awk  'NR==FNR{a[$2]++;next} !($2 in a)' ${maf_path}/MAF_5percent.txt ${input_path}/o${file}.out.fastGWA > o${file}_${analysis}_filter_5pMAF.txt

# Step 2: Now extract min p-value from the above filtered files 
# i) extract MAF 0.1%
awk '$10!=""{print $0}' o${file}_${analysis}_filter_01pMAF.txt | awk 'NR==1 {min=$10; row=$0} $10<min {min=$10; row=$0} END {print row}' > out_p${file}_${analysis}_min_pval_01percent.addGWA
# ii) extract MAF 0.5%
awk '$10!=""{print $0}' o${file}_${analysis}_filter_05pMAF.txt | awk 'NR==1 {min=$10; row=$0} $10<min {min=$10; row=$0} END {print row}' > out_p${file}_${analysis}_min_pval_05percent.addGWA
# iii) extract min p-value MAF 1%
awk '$10!=""{print $0}' o${file}_${analysis}_filter_1pMAF.txt  | awk 'NR==1 {min=$10; row=$0} $10<min {min=$10; row=$0} END {print row}' > out_p${file}_${analysis}_min_pval_1percent.addGWA
# iv) extract min p-value MAF 5%
awk '$10!=""{print $0}' o${file}_${analysis}_filter_5pMAF.txt  | awk 'NR==1 {min=$10; row=$0} $10<min {min=$10; row=$0} END {print row}' > out_p${file}_${analysis}_min_pval_5percent.addGWA

# cp the output
rsync -av *.addGWA $output_path

cd ..
#rm -rf $workdir
