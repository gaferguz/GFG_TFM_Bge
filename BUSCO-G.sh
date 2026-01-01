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

#busco -m genome -i FCS_GX_Cleaned_Bge2.fasta -o Busco_Bge2 -l insecta_odb12
#busco -m genome -i HCov.scaffolds.FINAL.fasta -o Busco_Sdor -l insecta_odb12
#busco -m genome -i NR_95_rep_seq.okay.mrna -o Busco_Transcriptome -l insecta_odb12
#busco -m genome -i TSA_c50_t1_v2_AMPIMD.fas -o Busco_Silva -l insecta_odb12
#busco -m genome -i RNASpades_tx.fasta -o Busco_Spa -l insecta_odb12
#busco -m genome -i NR_95_perlib_TPM_rep_seq.okay.mrna -o Busco_perLIB_tmp_25 -l insecta_odb12
#busco -m genome -i TPM1_10c_mrna.fasta -o Busco_TPM1_10c -l insecta_odb12
#busco -m genome -i TPM1_20c_mrna.fasta -o Busco_TPM1_20c -l insecta_odb12
#busco -m genome -i TPM1_subgr_mrna.fasta -o Busco_TPM1_subg -l insecta_odb12
#busco -m genome -i TPM1_avg_mrna.fasta -o Busco_TPM1_avg -l insecta_odb12
busco -m genome -i AMPs_Enriched_Final_TPM1_okaycull.mrna -o Busco_TPM1_enriched -l insecta_odb12
