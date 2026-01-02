#!/bin/bash
#SBATCH --job-name=Bge_Trinity
#SBATCH --qos=long-mem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=395GB
#SBATCH --time=3-00:00:00
#SBATCH --output=Bge_Trinity_10.11.2025.log


Trinity --seqType fq --samples_file samples.tsv --SS_lib_type RF --max_memory 300G --CPU 12
