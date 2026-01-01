#!/bin/bash
#SBATCH --job-name=trim_nov21_Nextera
#SBATCH --qos=long-mem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100GB
#SBATCH --time=5-00:00:00
#SBATCH --output=trim_nov21_Nextera_09.11.25.log

#!/bin/bash

echo ""
echo "------------------------------------------"
echo " SCRIPT ACTUALIZADO DE FILTRADO DE READS  "
echo "------------------------------------------"
echo ""
echo "Este scritp consta de tres pasos efectivos:"
echo ""
echo "--------------------------------------------------------------------"
echo "| 1. Decontaminacion de reads pertenecientes a Humano y a PhiX      |"
echo "| 2. Filtrado laxo con Trimmomatic para adaptadores Nextera         |"
echo "| 3. Filtrado con TrimGalore optimizado para adaptadores NextSeq550 |"
echo "---------------------------------------------------------------------"
echo ""

path=$(pwd)
format="fastq.gz"
resp="s" # pareadas? (s/n)


mkdir ${path}/Rcorrected
mkdir ${path}/TM_reads

if [ ${resp} == "s" ]
then
	echo -e "\tRaw_R1\tRaw_R2\tRaw_R1_nt\tRaw_R2_nt\tR1_Trimmomatic\tR2_Trimmommatic\tR1_Trimmomatic_nt\tR2_Trimmommatic_nt\tR1_Rcorrector\tR2_Rcorrector\tR1_Rcorrector_nt\tR2_Rcorrector_nt" | tr "\t" "\n" > Saneado_Results.txt
	for i in $(ls -1 ${path}/*.${format} | sed "s/.*\///g" | sed "s/..${format}//g" | sort | uniq)
	do
		read1=$(ls -1 ${path}/${i}* | grep "1.${format}" | sed "s/.*\///g")
		read2=$(ls -1 ${path}/${i}* | grep "2.${format}" | sed "s/.*\///g")
		echo "--------------------------------------------------------------------------------------"
		echo "Saneando: ${read1} + ${read2}"
		# Conteos del archivo crudo
		Raw_R1=$(zcat -f ${path}/${read1} | grep "^@" | wc -l)
		Raw_R2=$(zcat -f ${path}/${read2} | grep "^@" | wc -l)
		Raw_R1_nt=$(zcat -f ${path}/${read1} | awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}')
		Raw_R2_nt=$(zcat -f ${path}/${read2} | awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}')
		# Filtrado con Trimmomatic basico y poco restrictivo
		echo ""
		echo "Filtrado con Trimmomatic"
		trimmomatic PE -threads 8 ${path}/${read1} ${path}/${read2} \
			       ${path}/TM_reads/${i}1_trD_1.fastq ${path}/TM_reads/${i}1_un.trimmed.fastq \
			       ${path}/TM_reads/${i}2_trD_2.fastq ${path}/TM_reads/${i}2_un.trimmed.fastq \
			       ILLUMINACLIP:/home/gaferguz/Extras/adaptors/NexteraPE-PE.fa:2:30:10:8:True  SLIDINGWINDOW:4:20 MINLEN:50 LEADING:28
		# Conteos del archivo filtrado con Trimmomatic
		R1_Trimmomatic=$(zcat -f ${path}/TM_reads/${i}1_trD_1.fastq | grep "^@" | wc -l)
		R2_Trimmomatic=$(zcat -f ${path}/TM_reads/${i}2_trD_2.fastq | grep "^@" | wc -l)
		R1_Trimmomatic_nt=$(zcat -f ${path}/TM_reads/${i}1_trD_1.fastq | awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}')
		R2_Trimmomatic_nt=$(zcat -f ${path}/TM_reads/${i}2_trD_2.fastq | awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}')
		# Filtrado exhaustivo con trimgalore
		echo ""
		echo "Uso de Rcorrector"
		run_rcorrector.pl -1 ${path}/TM_reads/${i}1_trD_1.fastq -2 ${path}/TM_reads/${i}2_trD_2.fastq -od ./Rcorrected -t 12 -k '23' -maxcorK '4' -wk '0.95' -ek '100000000'
	 	# Conteos del archivo filtrado con Trimmomatic
		R1_Rcorrector=$(zcat -f ${path}/Rcorrected/${i}1_trD_1.cor.fq | grep "^@" | wc -l)
		R2_Rcorrector=$(zcat -f ${path}/Rcorrected/${i}2_trD_2.cor.fq | grep "^@" | wc -l)
		R1_Rcorrector_nt=$(zcat -f ${path}/Rcorrected/${i}1_trD_1.cor.fq | awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}')
		R2_Rcorrector_nt=$(zcat -f ${path}/Rcorrected/${i}2_trD_2.cor.fq | awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}')

	 	# Adicion de nueva linea al log de saneamiento
	 	echo -e "${i}\t${Raw_R1}\t${Raw_R2}\t${Raw_R1_nt}\t${Raw_R2_nt}\t${R1_Trimmomatic}\t${R2_Trimmomatic}\t${R1_Trimmomatic_nt}\t${R2_Trimmomatic_nt}\t${R1_Rcorrector}\t${R2_Rcorrector}\t${R1_Rcorrector_nt}\t${R2_Rcorrector_nt}" | tr "\t" "\n" > new_column.txt
		awk 'BEGIN {OFS="\t"} NR==FNR {a[FNR]=$1; next} {print $0, a[FNR]}' new_column.txt Saneado_Results.txt > Saneado_Results_${i}.tmp
		mv Saneado_Results_${i}.tmp Saneado_Results.txt
	done
	cd ${path} # Para no crear archivos tar con la ruta completa de directorios
fi




