#!/bin/bash
#SBATCH --job-name=Busco-TPM1_avg
#SBATCH --qos=long-mem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=300GB
#SBATCH --time=2-00:00:00
#SBATCH --output=Busco_TPM1_avg.log

# module load biotools
# module load anaconda
# conda activate busco

busco -m genome -i AMPs_Enriched_Final_TPM1_okaycull.mrna -o Busco_TPM1_enriched -l insecta_odb12
