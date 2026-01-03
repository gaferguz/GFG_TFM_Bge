# GFG_TFM_Bge
This repository is meant to store all the scripts used for **De novo transcriptome assembly**, Diffential Expression Analysis (**DET**), Weigthed Gene Correlation Network Analysis (**WGCNA**) and Gene Set Enrichment Analysis (**GSEA**) in *Blatella germanica* samples as part of the Master's Degree Final Project. All code is based on bash and R.

** **

Este repositorio reúne todos los scripts generados en el ensamblaje de un transcriptoma de novo, el análisis de expresión diferencial, análisis de redes y análisis de enriquecimiento funcional en *Blatella germanica* como parte del trabajo de fin de master.


# Pipeline:
```
- Quality control of fastq (Trimming/fasQC.sh)

- Assembly Trinity + RNASpades (Trimming/Sanear_[X]_TFM.sh)

- Pseudo-alignment using Kallisto (Alinear_kallisto.sh)

- TPM >= 1 filtering (Filtrado_TPM/TPM_Filtering.R)

- Creating a consensus (EvidentialGene/Generar_Consenso.sh)

- Merge okay and cull set (Filtrado_TPM/Seleccion_transcritos_okaycull.R)

- Pseudo-aligment using Kallisto (Alinear_kallisto.sh)

- Filtering based on TPM >=1 (Filtrado_TPM/TPM_Filtering.R)

- Final round of Kallisto on final transcriptome (Alinear_kallisto.sh)

- DEG analysis using DESeq2 + GSEA using clusterProfiler and Eggnog (Analisis_DET_GSEA.R)

- Weighted gene correlation network analysis (WGCNA/WGCNA_TFM.R)
```
