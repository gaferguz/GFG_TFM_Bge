#!/bin/bash
#SBATCH --job-name=Kallisto
#SBATCH --qos=long-mem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=20GB
#SBATCH --time=4-00:00:00
#SBATCH --output=Bge_Kallisto_Enriched_18.11.log

mkdir Kallisto # Creamos carpeta donde volcar todos los datos

# Definimos las variables principales
#index='/home/gaferguz/Blatella/Fused_9.11.2025/Trinity_idx'
#index='/home/gaferguz/Blatella/Fused_9.11.2025/SPADES_idx'
#index='/home/gaferguz/Blatella/Fused_9.11.2025/okaycull_idx'
index='/home/gaferguz/Blatella/Fused_9.11.2025/Enriched'
# Asumimos Pair_end acabados en R1/R2.${format}
format='cor.fq'

loc=$(pwd)

for i in $(ls *${format} | grep "R1" | sed "s/R1.*//g")
do
	Pair1=$(ls -1 ${i}* | grep "R1")
	Pair2=$(ls -1 ${i}* | grep "R2")
	outname=$(echo "${i}" | sed "s/\\.//g")
	echo "kallisto on samples ${i}"
	kallisto quant -i ${index} -o ${loc}/output_${outname} ${Pair1} ${Pair2} --rf-stranded -t 12  # --plaintext -t 12 # Quitamos el plaintext
	if [ ! -f ${loc}/Kallisto_Merged/Abundance.tsv ]

	then
    		cat ${loc}/output_${outname}/abundance.tsv | awk '{print $1,$4}' | sed "s/est_counts/est_counts_${outname}/g" >> Kallisto_Merged/Abundance_counts.tsv
    		cat ${loc}/output_${outname}/abundance.tsv | awk '{print $1,$NF}' | sed "s/tpm/tpm_${outname}/g" >> Kallisto_Merged/Abundance.tsv
	else
        	cat ${loc}/output_${outname}/abundance.tsv | awk '{print $NF}' | sed "s/tpm/tpm_${outname}/g" > temp.txt
        	awk 'NR==FNR {a[NR]=$NF; next} {print $0, a[FNR]}' temp.txt Kallisto_Merged/Abundance.tsv > temp2.txt
		mv temp2.txt Kallisto_Merged/Abundance.tsv

            	cat ${loc}/output_${outname}/abundance.tsv | awk '{print $4}' | sed "s/est_counts/est_counts_${outname}/g" > temp3.txt
                awk 'NR==FNR {a[NR]=$NF; next} {print $0, a[FNR]}' temp3.txt Kallisto_Merged/Abundance_counts.tsv > temp4.txt
                mv temp4.txt Kallisto_Merged/Abundance_counts.tsv


        fi
	echo "-------------------------------------------------------------------------------"
done

rm temp.txt
rm temp3.txt
