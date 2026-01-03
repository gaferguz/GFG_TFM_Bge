# Cargamos las librerias necesarias
library(dplyr)
library(ggplot2)
library(data.table)
library(MatrixGenerics)
library(matrixStats)
library(tximport)
library(DESeq2)
library(btools)
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
# remotes::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)
library(phyloseq)
# devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
library(ggpubr)
library(tidyverse)
library(clusterProfiler)


##################################
#           FUNCIONES            #
##################################


Deseq2_Vulcano <- function(deseq2_table, deseq2_plot, plot_outname, subtitle){
  
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = deseq2_table
  title_plot = deseq2_plot
  title_plot = str_replace(title_plot, "DESeq2_", "")
  title_plot = str_replace(title_plot, "_vs_", " vs ")
  title_plot = str_replace(title_plot, "_plot", "")
  
  plot = EnhancedVolcano(report,
                         title = title_plot,
                         lab = rownames(report),
                         selectLab = FALSE,
                         titleLabSize = 15,
                         subtitle = subtitle,
                         subtitleLabSize = 15,
                         captionLabSize = 15,
                         x = "log2FoldChange",
                         y = "padj",
                         FCcutoff = 1,                          
                         labSize = 7,
                         # legend = NA,
                         boxedLabels = TRUE,
                         pointSize = 2.0,
                         # axisLabSize = 10,
                         max.overlaps = 5,
                         legendPosition = "none",
                         gridlines.major = TRUE,
                         gridlines.minor = TRUE,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         pCutoff = 0.05) 
  
  ggsave(plot = plot, filename = plot_outname, height = 5, width = 4, dpi = 300)
    
}
 


Functional_Enrichment <- function(deseq2_table, deseq2_name, outdir, outdir_imm){
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = deseq2_table
  title = deseq2_name
  title = str_replace(title, "DESeq2_", "")
  title = str_replace(title, "_plot", "")
  
  reportf = unique(report)
  
  ############################################################################
  # GSEA
  # This time I am not filtering the changes
  gsea_list <- reportf %>%
    # dplyr::filter(log2FoldChange >= 1 & padj <= 0.05) %>%
    dplyr::arrange(desc(log2FoldChange))
  
  gsea_input <- gsea_list %>%
    dplyr::select(log2FoldChange) %>%
    unlist() %>%
    as.vector()
  names(gsea_input) <- rownames(gsea_list)
  
  # Gene enrichment analysis using the whole dataset sorted by FC
  enrichment_gsea <- GSEA(geneList = gsea_input,
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name,
                          minGSSize = 10,
                          maxGSSize = 500,
                          eps = 1e-10,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
  #save the enrichment result
  write.csv(file = paste0(outdir, "/", title, "_GSEA.csv"),
            x = enrichment_gsea@result)
  
  if (any(enrichment_gsea@result$p.adjust <= 0.05)){
    p <- ridgeplot(enrichment_gsea,
                   core_enrichment= TRUE, # Cambiado, en FALSE no existe
                   fill="p.adjust",
                   orderBy = "NES",
                   showCategory=100) +
      ggtitle("Ridge plot for GSEA")
    
    ggsave(filename = paste0(outdir, "/", title, "_GSEA_ridgeplot.pdf"),
           plot =  p,  dpi = 300, width = 21, height = 70, units = "cm")
  }
  # Filtering
  enrichment_gsea_immuno <- enrichment_gsea@result
  enrichment_gsea_immuno <- enrichment_gsea_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd",
                                                         enrichment_gsea_immuno$Description,
                                                         ignore.case = TRUE),]
  
  enrichment_gsea_immuno <- enrichment_gsea_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell|recombination",
                                                          enrichment_gsea_immuno$Description,
                                                          ignore.case = TRUE),]
  
  enrichment_gsea_immuno <- enrichment_gsea_immuno %>% filter(p.adjust < 0.05)
  
  if (nrow(enrichment_gsea_immuno) > 0) {
    write.csv(file = paste0(outdir_imm, "/", title,"_Immuno_GSEA.csv"), x = enrichment_gsea_immuno)
  }
} 
  

DifExp_pairwise <- function(mifactor, muestras, objeto_dds){
  
  # Inicializamos el dataframe para las listas de respuesta humoral
  Inmuno <- data.frame(matrix(ncol = 9, nrow = 0))
  x <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Group","TxID")
  colnames(Inmuno) <- x
  
  # Creamos una carpetas de salida
  outpath <- paste0(base::getwd(), "/", mifactor, "/")
  base::dir.create(outpath)
  outpath_vulcano <- paste0(outpath, "Vulcano_plots")
  base::dir.create(outpath_vulcano)
  outpath_tables <- paste0(outpath, "DESeq2_tables")
  base::dir.create(outpath_tables)
  outpath_GSEA <- paste0(outpath, "GSEA")
  base::dir.create(outpath_GSEA)
  outpath_GSEA_Immuno <- paste0(outpath, "GSEA_Immuno")
  base::dir.create(outpath_GSEA_Immuno)
  
  # Establecemos los grupos y las comparciones en funcion del factor
  lev <- levels(as.factor(muestras[,mifactor])) # get the variables
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
  
  for (j in 1:length(L.pairs)) {

    # Extraccion de los 2 factores de cada pareja posible
    factor1 = L.pairs[[j]][1]
    factor2 = L.pairs[[j]][2]
    
    # Creacion del nombre de la comparcion e inicalizacion del dataframe de resultados
    nam <- paste0("DESeq2_",factor1,"_vs_",factor2)
    plot_n <- paste0(nam,"_plot")
    dds_df <- data.frame(matrix(ncol = 7, nrow = 0))
    x <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","rank")
    colnames(dds_df) <- x
    
    # Aplicacion del test de Wald de DESeq2 con correccion "ash" del log2FoldChange
    factor1vs2 <- results(objeto_dds, contrast=c(mifactor, factor1, factor2), alpha = 0.05)
    factor1vs2 <- lfcShrink(objeto_dds, contrast=c(mifactor, factor1, factor2), res=factor1vs2, type = "ashr") #type = "normal") # Este tipo funciona
    factor1vs2 = data.frame(factor1vs2)
    dds_df = rbind(dds_df, factor1vs2)
    dds_df <- data.frame(dds_df)
    plot_n <- as.character(plot_n)
    
    # Filtrado por significancia y por listas IMD/Toll/AMPs
    reportf = unique(dds_df)
    reportf = reportf %>% filter(padj < 0.05)
    reportf = reportf %>% filter(abs(log2FoldChange) >= 1)
    subtitle_plot <- paste0("DEGs: ", nrow(reportf))
    selected_tx <- c()
    selected_tx <- c(selected_tx, rownames(reportf[grepl("_", rownames(reportf)),]))
    report_Inmuno <- reportf[grepl("_", rownames(reportf)),] # Usamos el patron "_"
    # report_Inmuno$Group <- rep(title_plot, nrow(report_Inmuno))
    
    # Usaremos este filtrado de lo que nos interesa para señalar los DEGs inmunos en el volcano plot
    report_Inmuno = unique(report_Inmuno)
    report_Inmuno$TxID <- rownames(report_Inmuno)
    Inmuno <- rbind(Inmuno, report_Inmuno)
    
    # Creamos un Vulcano Plot 
    Vulcname <- paste0(plot_n, ".png")
    Deseq2_Vulcano(dds_df, plot_n, paste0(outpath_vulcano,"/",Vulcname), subtitle_plot)
    
    # Creamos la tabla DESeq2 con resultados
    write.table(dds_df, file = paste0(outpath_tables,"/", gsub("_plot","_table", plot_n)), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Analisis GSEA con la funcion definida
    Functional_Enrichment(dds_df, plot_n, outpath_GSEA, outpath_GSEA_Immuno)
    
  }
  
  # Creamos el archivo con los datos de las listas IMD/Toll/AMPs
  write.table(Inmuno, file = paste0(outpath,"/Inmuno_DESeq2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}



##################################
#              MAIN              #
##################################

# Cargamos los archivos de anotacion funcional
term2gene <- read_tsv("./EggNOG/term2gene_GO.tsv")
term2name <- read_tsv("./EggNOG/term2name_GO.tsv")


# Listamos los archivos y generamos la ruta usando el patron output
# IMPORTANTE, LAS MUESTRAS ESTAN EN LA MATRIZ SEGUN SU ORDEN ALFABETICO EN LA VARIABLE
path <- "./Kallisto/"
x <- list.files(path)
x <- x[grepl("output", x)]
files <- file.path(path, x, "abundance.h5")
samples <- read.table("./samples.txt", header=TRUE)
rownames(samples) <- samples$sample

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

# # PARTE PARA AGREGAR POR GENES, NO SE HACE EN ESTE ANALISIS
# tx2gene <- data.frame(names(txi.kallisto$abundance[,1]), gsub("t[0-9].*","",(names(txi.kallisto$abundance[,1]))))
# colnames(tx2gene) <- c("txID","geneID")
# tx2gene$geneID <- gsub("_sm","_smt3", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
# txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ group) 
dds <- DESeq(dds)
vsd <- vst(dds)
PCAtissue <- plotPCA3D(vsd, intgroup = "group")

## Iteramos la funcion que hemos creado para los dos niveles de agregacion

factores <- c("tissue","condition")

for (i in factores) {
  des <- as.formula(paste0("~ ", i))
  dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = des) 
  dds <- DESeq(dds)
  # Si es por la columna de condicion, hacemos un objeto deseq2 por cada tejido
  if (i == "condition") {
    for (j in unique(samples$tissue)) {
      subsamples <- samples[grepl(j, samples$tissue),]
      sf <- paste(subsamples$sample, collapse="|")
      subfiles <- files[grepl(sf, files)]
      sub_txi.kallisto <- tximport(subfiles, type = "kallisto", txOut = TRUE)
      dds <- DESeqDataSetFromTximport(sub_txi.kallisto, colData = subsamples, design = des) 
      dds <- DESeq(dds)
      DifExp_pairwise(i, subsamples, dds)
    }
  }else{
    # En caso contratio, hacemos un test deseq2 en el que comparamos todos los tejidos
    DifExp_pairwise(i, samples, dds)
  }
}


