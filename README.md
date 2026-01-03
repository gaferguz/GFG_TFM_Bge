# GFG_TFM_Bge
This repository is meant to store all the scripts used for **De novo transcriptome assembly**, Diffential Expression Analysis (**DET**), Weigthed Gene Correlation Network Analysis (**WGCNA**) and Gene Set Enrichment Analysis (**GSEA**) in *Blatella germanica* samples as part of the Master's Degree Final Project. All code is based on bash and R.

** **

Este repositorio reúne todos los scripts generados en el ensamblaje de un transcriptoma de novo, el análisis de expresión diferencial, análisis de redes y análisis de enriquecimiento funcional en *Blatella germanica* como parte del trabajo de fin de master.


# Pipeline:

1 - Control de calidad de los .fastq (Trimming/fasQC.sh)

2 - Ensamblaje Trinity + RNASpades (Trimming/Sanear_[X]_TFM.sh)

3 - Pseudo-alineamiento usando Kallisto (Alinear_kallisto.sh)

4 - Filtro por TPM >= 1 (Filtrado_TPM/TPM_Filtering.R)

5 - Creando un consenso (EvidentialGene/Generar_Consenso.sh)

6 - Fusion de transcritos "okay" y "cull" (Filtrado_TPM/Seleccion_transcritos_okaycull.R)

7 - Pseudo-alineamiento usando Kallisto (Alinear_kallisto.sh)

8 - Filtro por TPM >= 1 (Filtrado_TPM/TPM_Filtering.R)

9 - Ronda final de Kallisto en el transcriptoma definitivo (Alinear_kallisto.sh)

10 - Analysis de expresion diferencial con DESeq2 + GSEA con clusterProfiler y Eggnog (Analisis_DET_GSEA.R)

11 - Weighted gene correlation network analysis (WGCNA/WGCNA_TFM.R)

El script para procesar las anotaciones del transcriptoma están en la carpeta /EggNog
