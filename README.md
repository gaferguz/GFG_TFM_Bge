# GFG_TFM_Bge
This repository is meant to collect a number of scripts used for De novo transcriptome assembly, diffential expression analysis, Weigthed Gene Correlation Network Analysis (WGCNA) and Gene Set Enrichment Analysis (GSEA) as part of the final project of a masters degree thesis. 

# Scripts of the pipeline (bash and R scripting):
- Trimming and sanitizing raw fastq files with Trimmomatic + Rcorrector
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



