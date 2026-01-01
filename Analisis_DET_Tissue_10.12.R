library(dplyr)
library(ggplot2)
library(data.table)
library(MatrixGenerics)
library(matrixStats)
library(tximport)
library(DESeq2)
library(btools)

###############

# Parte copiada del TFM de Nuria Andres / bioconductor vignetes:
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

# Listamos los archivos y generamos el path, en este caso la tabla se ha hecho en excel
path <- "./Kallisto_Enriched/"
x <- list.files(path)
x <- x[grepl("output", x)]
# Es un outlier que habiamos visto pero al agregar por gene ya no lo veo mas
# x <- x[!grepl("Control-Malpighi-3", x)] 


files <- file.path(path, x, "abundance.h5")

# IMPORTANTE LAS SAMPLES ESTAN PUESTAS SEGUN VAN ORDENADAS ALFABEITCAMENTE EN LA VARIABLE
# x <- list.files(path), antes lo hacia segun el orden de ls en linux y las columnas estaban
# mal anotadas, REVISAR SIEMPRE BIEN SI LOS COUNTS COINCIDEN con las abundancias individuales
# de cada archivo
samples <- read.table("./samples.txt", header=TRUE)
rownames(samples) <- samples$sample

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE) # Esta da error

# # PARTE DE DEGS QUE AQUI SE EVITA (01.12.2025)
# tx2gene <- data.frame(names(txi.kallisto$abundance[,1]), gsub("t[0-9].*","",(names(txi.kallisto$abundance[,1]))))
# colnames(tx2gene) <- c("txID","geneID")
# tx2gene$geneID <- gsub("_sm","_smt3", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
# txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
# 
# # Original
# # txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ tissue) 
dds <- DESeq(dds)
vsd <- vst(dds)
PCAtissue <- plotPCA3D(vsd, intgroup = "tissue")


# Ajuste al modelo de dispersion (test de Ines)
# dds_test <- collapseReplicates(dds, dds$condition)
# plotDispEsts(dds_test, main="Dispersion plot")
# 
# min(PCAcond$pca$rotation[,"PC1"])
# PCAcond$pca$rotation[,"PC1"][grep("CG", names(PCAcond$pca$rotation[,"PC1"]))]
# PCAcond$pca$rotation
# pca = prcomp(t(assay(vsd1)))
# PCAcond$pca$x[, 1]
# assay(vsd1)
# pca$x


# Vamos a iterar por tejidos y crear un nuevo objeto deseq por cada subconjunto
# argumentamos que de los contrario se inflan los valores de logFC en base a 
# una basemean y un factor de normalizacion con grupos muy distintos, sabiendo
# que las muestras se segregan bien por tejido y que hay suficientes replicas

deseq2_tables <- list()
deseq2_plots <- list()


contador = 0
lev2 <- levels(as.factor(samples[,2])) # get the variables
L.pairs <- combn(seq_along(lev2), 2, simplify = FALSE, FUN = function(i) lev2[i])
# sublist <- L.pairs[grepl("_Ctrl", L.pairs)]
  
for (j in 1:length(L.pairs)) {
    contador = contador + 1
    factor1 = L.pairs[[j]][1]
    factor2 = L.pairs[[j]][2]
    nam <- paste0("DESeq2_",factor1,"_vs_",factor2)
    plot_n <- paste0(nam,"_plot")
    dds_df <- data.frame(matrix(ncol = 7, nrow = 0))
    x <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","rank")
    colnames(dds_df) <- x
    dds_all = dds
    
    factor1vs2 <- results(dds_all, contrast=c("tissue", factor1, factor2), alpha = 0.05)
    factor1vs2 <- lfcShrink(dds_all, contrast=c("tissue", factor1, factor2), res=factor1vs2, type = "ashr") #type = "normal") # Este tipo funciona
    # factor1vs2 <- lfcShrink(dds_all, coef=paste0("condition_", factor1, "_vs_",factor2), res=factor1vs2, type = "apeglm") 
    # factor1vs2 = factor1vs2 %>% data.frame()
    factor1vs2 = data.frame(factor1vs2)
    dds_df = rbind(dds_df, factor1vs2)
    assign(nam,dds_df)
    deseq2_tables[[contador]] <- dds_df
    deseq2_plots[[contador]] <- as.character(plot_n)
    # command <- paste0("echo ",nam," >> ./DE.txt")
    # system(command)
    # system("echo \n >> ./DE.txt")
}

# contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
# res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
# res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken)

# The order of the names determines the direction of fold change that is reported. 
# The name provided in the second element is the level that is used as baseline. 
# So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe relative to the control.

#####################
# VULCANOS - DESeq2 #
#####################
library(tidyverse)
library(phyloseq)
# devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
library(ggpubr)


deseq2_plots_objects <- list()
Inmnuno <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Group","TxID")
colnames(Inmnuno) <- x


for(i in 1:length(deseq2_plots)){ 
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = data.frame(deseq2_tables[[i]])
  title_plot = deseq2_plots[[i]]
  title_plot = str_replace(title_plot, "DESeq2_", "")
  title_plot = str_replace(title_plot, "_vs_", " vs ")
  title_plot = str_replace(title_plot, "_plot", "")
  
  # Filtrdo por significancia 
  reportf = unique(report)
  reportf = report %>% filter(padj < 0.05)
  # reportf = reportf %>% filter(abs(log2FoldChange) >= 2.3219) # Lo2FC referencia de Nuria Andres
  reportf = reportf %>% filter(abs(log2FoldChange) > 1) # Filtro de Silva et al 2024
  # reportf <- reportf[order(reportf$log2FoldChange, decreasing = TRUE),]
  # selected_tx <- head(rownames(reportf), n=5) # Top 5 por arriba y abajo
  # reportf <- reportf[order(reportf$log2FoldChange, decreasing = FALSE),]
  # selected_tx <- c(selected_tx, head(rownames(reportf), n=5))
  subtitle_plot <- paste0("DEGs: ", nrow(reportf))
  selected_tx <- c()
  # Ahora buscamos los que estan anotados a parte y los marcamos tambien
  # selected_tx <- c(selected_tx, rownames(reportf[!grepl("NonamEVm", rownames(reportf)),])) # Antigua
  selected_tx <- c(selected_tx, rownames(reportf[grepl("_", rownames(reportf)),]))
  # report_Inmuno <- report[!grepl("NonamEVm", rownames(report)),] # Antigua
  report_Inmuno <- reportf[grepl("_", rownames(reportf)),] # Ahora aprovechamos _
  report_Inmuno$Group <- rep(title_plot, nrow(report_Inmuno))
  
  # Usaremos este filtrado de lo que nos interesa para señalar los DEGs inmunos en el volcano plot
  report_Inmuno = unique(report_Inmuno)
  report_Inmuno$TxID <- rownames(report_Inmuno)
  # report_Inmuno = report_Inmuno %>% filter(padj < 0.05)
  # reportf = reportf %>% filter(abs(log2FoldChange) >= 2.3219) # Lo2FC referencia de Nuria Andres
  # report_Inmuno = report_Inmuno %>% filter(abs(log2FoldChange) >= 2)
  Inmnuno <- rbind(Inmnuno, report_Inmuno)
  
  plot = EnhancedVolcano(report,
                         title = title_plot,
                         # lab = NA,
                         lab = rownames(report),
                         selectLab = selected_tx,
                         titleLabSize = 20,
                         subtitle = subtitle_plot,
                         subtitleLabSize = 20,
                         # caption = NULL,
                         captionLabSize = 20,
                         x = "log2FoldChange",
                         y = "padj",
                         FCcutoff = 1, # LogFoldChange mayor a 2 (referencia de Silva et al 2024)
                         # FCcutoff = 2.3219, # LogFoldChange mayor a 5 (referencia de Nuria Andres)
                         labSize = 5,
                         # legend = NA,
                         boxedLabels = TRUE,
                         # parseLabels = TRUE,
                         # labCol = 'black',
                         # labFace = 'bold',
                         pointSize = 2.0,
                         # legendLabSize = 0.001,
                         # legendIconSize = 0.00001,
                         # axisLabSize = 10,
                         max.overlaps = 5,
                         legendPosition = "none",
                         gridlines.major = TRUE,
                         gridlines.minor = TRUE,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         pCutoff = 0.05) 
  
  # Assing all objects to its variable name
  assign(as.character(deseq2_plots[[i]]), plot)
  deseq2_plots_objects[[i]] <-  plot
}

write.table(Inmnuno, file = paste0(path,"/Inmuno_DESeq2_res.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# y <- ggarrange(plotlist = deseq2_plots_objects, ncol = 6, nrow = 6)
library(gridExtra)
library(grid)

# Se tiene que hacer por separado porque todo junto supera los limites de agrupamiento grafico con estos objetos y el grid o el ggarrange

T1 <- grid.arrange(DESeq2_FB_vs_FG_plot, 
                   DESeq2_FB_vs_MG_plot,
                   DESeq2_FB_vs_HG_plot,
                   nrow=1,
                   ncol=3)

T2 <- grid.arrange(DESeq2_FB_vs_SG_plot, 
                   DESeq2_FB_vs_MT_plot,
                   DESeq2_FG_vs_HG_plot,
                   nrow=1,
                   ncol=3)

T3 <- grid.arrange(DESeq2_FG_vs_MG_plot, 
                   DESeq2_FG_vs_MT_plot,
                   DESeq2_FG_vs_SG_plot,
                   nrow=1,
                   ncol=3)

T4 <- grid.arrange(DESeq2_HG_vs_MG_plot, 
                   DESeq2_HG_vs_MT_plot,
                   DESeq2_HG_vs_SG_plot,
                   nrow=1,
                   ncol=3)

T5 <- grid.arrange(DESeq2_MG_vs_MT_plot, 
                   DESeq2_MG_vs_SG_plot,
                   DESeq2_MT_vs_SG_plot,
                   nrow=1,
                   ncol=3)


##############################################################################


for (i in 1:length(deseq2_plots)) {
  x <- get(gsub("_plot", "", deseq2_plots[[i]]))
  saveRDS(x, file = paste0(path,"/", gsub("_plot", "", deseq2_plots[[i]]), ".rds"))
}

for (i in 1:length(deseq2_plots)) {
  x <- get(deseq2_plots[[i]])
  saveRDS(x, file = paste0(path,"/", deseq2_plots[[i]], ".rds"))
}


##############################################################################
##############################################################################

#                           FUNCTIONAL ENRICHMENT

##############################################################################
##############################################################################

## Enriquecimiento funcional usando eggnog

# install.packages(c("BiocManager", "devtools"))
# library("devtools")
# library("BiocManager")
# BiocManager::install("tidyverse")
# devtools::install_github("GuangchuangYu/GOSemSim")
# devtools::install_github("GuangchuangYu/clusterProfiler")


library(devtools)
library(tidyverse)
library(clusterProfiler)
library(ontologyIndex)

# prepare the term to gene table
eggNOG <- read_tsv("/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/Bge_Annotations.tsv") %>%
  dplyr::select(GOs, query) %>%
  dplyr::filter(GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  mutate(gene = gsub("\\..*", "", query)) %>%
  select(GOs, gene) %>%
  distinct() %>%
  drop_na()
colnames(eggNOG) <- c("term", "gene")

# prepare the term to name table
# Tengo que descargarme de manera manual la DB de go.obo
# https://geneontology.org/docs/download-ontology/
ontology <- get_ontology(file = "/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/go.obo",
                         propagate_relationships = "is_a",
                         extract_tags = "everything",
                         merge_equivalent_terms = TRUE)

# ontology <- get_ontology(file = "/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/goslim_drosophila.obo",
#                          propagate_relationships = "is_a",
#                          extract_tags = "everything",
#                          merge_equivalent_terms = TRUE)

eggNOG_term <- eggNOG %>%
  mutate(name = ontology$name[term]) %>%
  select(c(term, name)) %>%
  distinct() %>%
  drop_na() %>%
  filter(!grepl("obsolete", name))

eggNOG <- eggNOG %>%
  filter(term %in% eggNOG_term$term)

# save the results
write_tsv(x = eggNOG, file = "/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/term2gene_GO.tsv")
write_tsv(x = eggNOG_term, file = "/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/term2name_GO.tsv")


# perform ORA
term2gene <- read_tsv("/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/term2gene_GO.tsv")
term2name <- read_tsv("/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/term2name_GO.tsv")


# KEGG ORA
eggNOG_kegg <- read_tsv("/media/gfg017/SEAGATE_GFG/9-TFM-2026/7.EggNOG/Bge_Annotations.tsv") %>%
  dplyr::select(KEGG_ko, query) %>%
  dplyr::filter(KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  dplyr::mutate(gene = gsub("\\..*", "", query)) %>%
  dplyr::mutate(term = gsub("ko:", "", KEGG_ko)) %>%
  dplyr::select(term, gene) %>%
  distinct() %>%
  drop_na()



deseq2_plots_objects <- list()

for(i in 1:length(deseq2_plots)){ 
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = data.frame(deseq2_tables[[i]])
  title = deseq2_plots[[i]]
  title = str_replace(title, "DESeq2_", "")
  
  reportf = unique(report)
  
  ############################################################################
  # ORA
  # read the gene list of interest
  over_set <- reportf %>%
    dplyr::filter(log2FoldChange >= 1 & padj <= 0.05)
  
  under_set <- reportf %>%
    dplyr::filter(log2FoldChange <= -1 & padj <= 0.05)
  
  over_set <- rownames(over_set)
  under_set <- rownames(under_set)
  
  over_enrichment <- enricher(over_set,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name,
                         pvalueCutoff = 0.05,
                         # universe = background_genes,
                         qvalueCutoff = 0.05)
  
  over_enrichment_immuno <- over_enrichment@result
  over_enrichment_immuno <- over_enrichment_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd", 
                                                         over_enrichment_immuno$Description, 
                                                         ignore.case = TRUE),]
  over_enrichment_immuno <- over_enrichment_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell", 
                                                          over_enrichment_immuno$Description, 
                                                          ignore.case = TRUE),]
  over_enrichment_immuno <- over_enrichment_immuno %>% filter(p.adjust < 0.05)
  if (nrow(over_enrichment_immuno) > 0) {
    write.csv(file = paste0(title,"_Immuno_over_ORA.csv"), x = over_enrichment_immuno)
  }
  
  #save the enrichment result
  write.csv(file = paste0(title,"_over_ORA.csv"), x = over_enrichment@result)
  
  if (any(over_enrichment@result$p.adjust <= 0.05)){
    p <- dotplot(over_enrichment,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
      ggtitle("dotplot for GO ORA")
    
    ggsave(filename = paste0(title,"_over_ORA_dotplot.pdf"),            
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  }
  
  under_enrichment <- enricher(under_set,
                              TERM2GENE = term2gene,
                              TERM2NAME = term2name,
                              pvalueCutoff = 0.05,
                              # universe = background_genes,
                              qvalueCutoff = 0.05)
  # Filtering
  under_enrichment_immuno <- under_enrichment@result
  under_enrichment_immuno <- under_enrichment_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd", 
                                                           under_enrichment_immuno$Description, 
                                                           ignore.case = TRUE),]
  under_enrichment_immuno <- under_enrichment_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell", 
                                                            under_enrichment_immuno$Description, 
                                                            ignore.case = TRUE),]
  under_enrichment_immuno <- under_enrichment_immuno %>% filter(p.adjust < 0.05)
  if (nrow(under_enrichment_immuno) > 0) {
    write.csv(file = paste0(title,"_Immuno_under_ORA.csv"), x = under_enrichment_immuno)
  }
  
  #save the enrichment result
  write.csv(file = paste0(title,"under_ORA.csv"), x = under_enrichment@result)
  
  if (any(under_enrichment@result$p.adjust <= 0.05)){
    p <- dotplot(under_enrichment,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
      ggtitle("dotplot for GO ORA")
    
    ggsave(filename = paste0(title,"_under_ORA_dotplot.pdf"),            
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  }
  
  ############################################################################
  # KEGG
  # read the gene list of interest
  # create a list of kegg ortholog that I am interested in
  over_set_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene %in% over_set) %>%
    unlist() %>%
    as.vector()
  
  under_set_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene %in% under_set) %>%
    unlist() %>%
    as.vector()
  
  # create a list of kegg ortholog that includes all kegg orthologs which form my background
  background_kegg <- eggNOG_kegg$term
  
  over_enrichment_kegg <- enrichKEGG(over_set_kegg,
                                organism = "ko",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                universe = background_kegg,
                                minGSSize = 10,
                                maxGSSize = 500,
                                qvalueCutoff = 0.05,
                                use_internal_data = FALSE)
  
  under_enrichment_kegg <- enrichKEGG(under_set_kegg,
                                organism = "ko",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                universe = background_kegg,
                                minGSSize = 10,
                                maxGSSize = 500,
                                qvalueCutoff = 0.05,
                                use_internal_data = FALSE)
  #save the enrichment result
  write.csv(file = paste0(title,"_over_enrichment_KEGG.csv"),                
            x = over_enrichment_kegg@result)
  
  write.csv(file = paste0(title,"_under_enrichment_KEGG.csv"),                
            x = under_enrichment_kegg@result)
  
  if (any(over_enrichment_kegg@result$p.adjust <= 0.05)){
    p <- dotplot(over_enrichment_kegg,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
      ggtitle("dotplot for KEGG ORA")
    
    ggsave(filename = paste0(title,"_over_enrichment_KEGG.pdf"),                # EDIT THIS
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  }
  
  if (any(under_enrichment_kegg@result$p.adjust <= 0.05)){
    p <- dotplot(under_enrichment_kegg,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
      ggtitle("dotplot for KEGG ORA")
    
    ggsave(filename = paste0(title,"_under_enrichment_KEGG.pdf"),                # EDIT THIS
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  }
  
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
  write.csv(file = paste0(title, "_GSEA.csv"),             
            x = enrichment_gsea@result)
  
  if (any(enrichment_gsea@result$p.adjust <= 0.05)){
    p <- ridgeplot(enrichment_gsea,
                   core_enrichment= TRUE, # Cambiado, en FALSE no existe
                   fill="p.adjust",
                   orderBy = "NES",
                   showCategory=100) +
      ggtitle("Ridge plot for GSEA")
    
    ggsave(filename = paste0(title, "_GSEA_ridgeplot.pdf"),                
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  }
  # Filtering
  enrichment_gsea_immuno <- enrichment_gsea@result
  enrichment_gsea_immuno <- enrichment_gsea_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd", 
                                                         enrichment_gsea_immuno$Description, 
                                                         ignore.case = TRUE),]
  
  enrichment_gsea_immuno <- enrichment_gsea_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell", 
                                                          enrichment_gsea_immuno$Description, 
                                                          ignore.case = TRUE),]
  
  enrichment_gsea_immuno <- enrichment_gsea_immuno %>% filter(p.adjust < 0.05)
  if (nrow(enrichment_gsea_immuno) > 0) {
    write.csv(file = paste0(title,"_Immuno_GSEA.csv"), x = enrichment_gsea_immuno)
  }
}















