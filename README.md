# GFG_TFM_Bge
This repository is meant to store all the scripts used for De novo transcriptome assembly, Diffential Expression Analysis, Weigthed Gene Correlation Network Analysis (WGCNA) and Gene Set Enrichment Analysis (GSEA) in **Blatella germanica** as part of the final project of a masters degree Final Project. All code is based on bash and R.

# Pipeline:
- Trimming and sanitizing raw fastq (Trimmomatic + Rcorrector)
- Quality control
- Assembly Trinity + RNASpades
- Pseudo-alignment using Kallisto
- TPM >= 1 filtering
- Creating a consensus and using Kallisto again
- Filtering based on TPM >=1
- Final round of Kallisto on final transcriptome
- DEG analysis on trasncripts using DESeq2
- Gene enrichment analysis on transcripts using gProfiler
- Weighted gene correlation network analysis



