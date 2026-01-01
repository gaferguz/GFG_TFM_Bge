#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --qos=long-mem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=100GB
#SBATCH --time=1-00:00:00
#SBATCH --output=FASTQC.log

fastqc -t 12 *.fq
