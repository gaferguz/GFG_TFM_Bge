# script to perform WGCNA
# setwd("~/Desktop/demo/WGCNA")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
# remotes::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)
library(gridExtra)
library(tximport)


# Listamos los archivos y generamos el path (igual que en el analisis DET + GSEA)
x <- list.files("./Kallisto/")
files <- file.path("./Kallisto/", x, "abundance.h5")
files <- files[grepl("output_", files)]
samples <- read.table(file.path("./samples.txt"), header=TRUE)
rownames(samples) <- samples$sample

# Creando clase usando counts de kallisto
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

# Objeto DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ condition) 
dds <- DESeq(dds)
data <- assay(dds)

# Normalizamos los datos usando la transformacion estabilizadora de la varianza
dds_norm <- vst(dds) 

# get normalized counts
norm.counts <- t(assay(dds_norm))

# Construccion de la red WGNA
# Creamos un vector para analizar distintos 
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Usamos la siguiente funcion para buscar el theshold que vamos aplicar la matriz de ayacencia
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)


sft.data <- sft$fitIndices

# Creamos un plot de visualizacion
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Poder', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Poder', y = 'Conectividad media') +
  theme_classic()
  
# Creamos un grafico para ver el valor al que se estabilizan ambos graficos (en mi caso 9)
grid.arrange(a1, a2, nrow = 2) 

# Convertimos los tipos de datos a una matriz numerica
norm.counts[] <- sapply(norm.counts, as.numeric)

# Establecemos los parametros
soft_power <- 9

# Funcion "all-in-one" que crea la matriz de adyacencia, la TOM y los modulos fusionados
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

# Sacamos los eigenGenes
module_eigengenes <- bwnet$MEs


# Asociamos los modulos a los factores de agrupacion de muestras

# Primero creamos columnas binarias para cada factor de agrupacion 
traits <- colData %>% 
  mutate(FB_state_bin = ifelse(grepl('FB', tissue), 1, 0)) %>% 
  select(5)

traits_tisgroup <- colData %>% 
  mutate(FBCtrl_state_bin = ifelse(grepl('FB_Ctrl', condition), 1, 0)) %>% 
  select(5)

# Binarizamos las variables categoricas
colData$tissue <- factor(colData$tissue, levels = c("FB","FG","HG","MT","MG","SG"))
colData$condition <- factor(colData$condition, levels = c("FB_Ctrl","FG_Ctrl","HG_Ctrl",
                                                          "MT_Ctrl","MG_Ctrl","SG_Ctrl",
                                                          "FB_G2r","HG_GFaFeces","SG_GFaFeces",
                                                          "HG_GFa","SG_GFa","FB_GFa",
                                                          "FG_GFa","FG_GFaFeces","MT_GFaFeces",
                                                          "MG_GFaFeces","MT_GFa","MG_GFa"))


tissue.out <- binarizeCategoricalColumns(colData$tissue,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)

condition.out <- binarizeCategoricalColumns(colData$condition,
                                        includePairwise = FALSE,
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)

traits <- cbind(traits, tissue.out)
traits_tisgroup <- cbind(traits_tisgroup, condition.out)


# Definimos el numero de muestras y genes
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Hacemos una correlacion de Pearson para los dos tipos de agrupacion (tejido y condicion+tejido)
module.trait.corr <- WGCNA::cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

module.traits_tg.corr <- WGCNA::cor(module_eigengenes, traits_tisgroup, use = 'p')
module.traits_tg.corr.pvals <- corPvalueStudent(module.traits_tg.corr, nSamples)

# Creamos heatmaps especificos

# TEJIDOS
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
heatmap.data <- heatmap.data %>% column_to_rownames(var = 'Row.names')

# Ajustamos manualmente los indices para las columnas que son del eje x e y
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[39:44],
             y = names(heatmap.data)[1:38],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)


# TEJIDOS + GRUPOS
heatmap.data.tg <- merge(module_eigengenes, traits_tisgroup, by = 'row.names')
heatmap.data.tg <- heatmap.data.tg %>% column_to_rownames(var = 'Row.names')

# Ajustamos manualmente los indices para las columnas que son del eje x e y
CorLevelPlot(heatmap.data.tg,
             x = names(heatmap.data.tg)[39:56],
             y = names(heatmap.data.tg)[1:38],
             col = c("blue1", "skyblue", "white", "pink", "red"))





# ANALISIS INTRAMODULO


# Calculamos MM (membresia al modulo) y p-valor asociado de cada transcrito al modulo

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
MMcors <- as.data.frame(module.membership.measure)

# Calculamos la significancia de cada transcrito en cada modulo y filtramos por significancia > 0.6
gene.signf.corr <- cor(norm.counts, traits$data.MT.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr_filt <- as.data.frame(gene.signf.corr)
gene.signf.corr_filt <-  as.character(rownames(gene.signf.corr_filt %>% dplyr::filter(.[[1]] > 0.6)))

# MODIFICACION DEL 21.12.2025
modNames = substring(names(bwnet$MEs), 3) #extract module names

# Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(norm.counts, bwnet$MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(norm.counts)))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(norm.counts, traits$FB_state_bin, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(norm.counts)))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")
head(GSPvalue)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(100)

# Choose Top Hub in each module function

colores <- bwnet$colors
Hubs <- chooseTopHubInEachModule(
  norm.counts, 
  colorh = colores, 
  # omitColors = "grey", 
  power = 9,
  type = "signed")

TopH <- topHubs(norm.counts, colorh = colores, power = 9, omitColors = NA)

gene.signf.corr_filt <- as.data.frame(gene.signf.corr)
gene.signf.corr_filt <-  as.character(rownames(gene.signf.corr_filt %>% dplyr::filter(.[[1]] > 0.6)))

# RASTREO A NIVEL DE TEJIDOS (21.12.2025)
# Ejemplo, se evita porque se pierda absolutamente todo
MM_MEred_SG <- MMcors[grepl("^MEred$", rownames(MMcors)),]
MM_MEred_SG <- as.data.frame(t(MM_MEred_SG))
MM_MEred_SG <- MM_MEred_SG %>% dplyr::filter(MEred > 0.83)
MM_MEred_SG <- as.character(rownames(MM_MEred_SG))
MM_MEred_SG <- MM_MEred_SG[MM_MEred_SG %in% gene.signf.corr_filt]

Hubs_MEred <- TopH[grepl("^red$", TopH$module),]
Hubs_MEred <- Hubs_MEred[order(Hubs_MEred$connectivity_rowSums_adj, decreasing = TRUE),]

paste(Hubs_MEred$gene[1:150], collapse = "|")

# Este modulo es enorme pero no muestra nada relacionado exclusivamente con respuesta inmune
Hubs_MEbrown_MTpos <- TopH[grepl("^brown$", TopH$module),]
Hubs_MEbrown_MTpos <- Hubs_MEbrown_MTpos[order(Hubs_MEbrown_MTpos$connectivity_rowSums_adj, decreasing = TRUE),]

paste(Hubs_MEbrown_MTpos$gene, collapse = "|")
paste(Hubs_MEbrown_MTpos$gene[1:150], collapse = "|")

# Vamos al otro que no tiene defensinas en el pero tiene una Tollo y es pequeño
Hubs_MEs3_MTpos <- TopH[grepl("^skyblue3$", TopH$module),]
Hubs_MEs3_MTpos <- Hubs_MEs3_MTpos[order(Hubs_MEs3_MTpos$connectivity_rowSums_adj, decreasing = TRUE),]

paste(Hubs_MEs3_MTpos$gene, collapse = "|")


# RASTREO A NIVEL DE GRUPOS 

# De esta que es la mas grande no ha salido nada
Hubs_MEma_FBGF_pos <- TopH[grepl("^magenta$", TopH$module),]
Hubs_MEma_FBGF_pos <- Hubs_MEma_FBGF_pos[order(Hubs_MEma_FBGF_pos$connectivity_rowSums_adj, decreasing = TRUE),]

paste(Hubs_MEma_FBGF_pos$gene, collapse = "|")

Hubs_MEmb_HGFeces_pos <- TopH[grepl("^midnightblue$", TopH$module),]
Hubs_MEmb_HGFeces_pos <- Hubs_MEmb_HGFeces_pos[order(Hubs_MEmb_HGFeces_pos$connectivity_rowSums_adj, decreasing = TRUE),]

paste(Hubs_MEmb_HGFeces_pos$gene, collapse = "|")



# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.


# the grey module is omitted
topHubs <- function (datExpr, colorh, omitColors = "grey", power = 2, type = "signed", 
                     ...) 
{
  # modified from chooseTopHubInEachModule, but return the table of all genes connectivity
  isIndex = FALSE
  modules = names(table(colorh))
  if (!is.na(omitColors)[1]) 
    modules = modules[!is.element(modules, omitColors)]
  if (is.null(colnames(datExpr))) {
    colnames(datExpr) = 1:dim(datExpr)[2]
    isIndex = TRUE
  }
  
  connectivity_table <- data.frame(matrix(ncol = 3)) %>% setNames(c('gene', 'connectivity_rowSums_adj', 'module'))
  hubs = rep(NA, length(modules))
  names(hubs) = modules
  for (m in modules) {
    adj = adjacency(datExpr[, colorh == m], power = power, 
                    type = type, ...)
    
    hub = which.max(rowSums(adj))
    
    hubs[m] = colnames(adj)[hub]
    
    sorted_genes <- rowSums(adj) %>% sort(decreasing = T) %>% as.data.frame()  %>%  
      tibble::rownames_to_column() %>% setNames(c('gene', 'connectivity_rowSums_adj')) %>% mutate(module = m)
    connectivity_table <- connectivity_table %>% rbind(sorted_genes)
    
    
    
  }
  if (isIndex) {
    hubs = as.numeric(hubs)
    names(hubs) = modules
  }
  return(connectivity_table %>% na.omit)
}



