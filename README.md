# GFG_TFM_Bge
This repository is meant to store all the scripts used for **De novo transcriptome assembly**, Diffential Expression Analysis (**DET**), Weigthed Gene Correlation Network Analysis (**WGCNA**) and Gene Set Enrichment Analysis (**GSEA**) in *Blatella germanica* samples as part of the Master's Degree Final Project. All code is based on bash and R.

Este repositorio reúne todos los scripts generados en el ensamblaje de un transcriptoma de novo, el análisis de expresión diferencial, análisis de redes y análisis de enriquecimiento funcional en *Blatella germanica* como parte del trabajo de fin de master.

** **

# Pipeline:
```
Trimming and sanitizing raw fastq (Trimmomatic + Rcorrector)
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
```
