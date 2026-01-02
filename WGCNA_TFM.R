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


allowWGCNAThreads()          # allow multi-threading (optional)


# Listamos los archivos y generamos el path, en este caso la tabla se ha hecho en excel
x <- list.files("../5.DEG/Kallisto_Enriched/")
files <- file.path("../5.DEG/Kallisto_Enriched/", x, "abundance.h5")
files <- files[grepl("output_", files)]
samples <- read.table(file.path("../5.DEG/samples.txt"), header=TRUE)
rownames(samples) <- samples$sample

# Creando clase usando counts de kallisto
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

# Objeto DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ condition) 
dds <- DESeq(dds)
data <- assay(dds)

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- WGCNA::goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2
data_orig <- data
data <- vst(data_orig)

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples (in this case we wont)
data.subset <- data_orig
# samples.to.be.excluded <- c('GFaFeces-Salivares-1', 'GFaFeces-HG-4')
# data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
# 

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- samples 
# colData <- phenoData %>% 
#   filter(!row.names(.) %in% samples.to.be.excluded)

#   
# # fixing column names in colData
# names(colData)
# names(colData) <- gsub(':ch1', '', names(colData))
# names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset)) # TRUE
all(rownames(colData) == colnames(data.subset)) # TRUE


# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ condition) # not spcifying model



# perform variance stabilization
dds_norm <- vst(dds) # we get a warning of fitType = 'parametric' not capturing well the dispersion trend


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
  

grid.arrange(a1, a2, nrow = 2) # We choose 12 as a soft theshold, over 0.8 R2 and minimize mean conectivity


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 9
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module



# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(FB_state_bin = ifelse(grepl('FB', tissue), 1, 0)) %>% 
  select(5)

traits_group <- colData %>% 
  mutate(Ctrl_state_bin = ifelse(grepl('Ctrl', group), 1, 0)) %>% 
  select(5)

traits_tisgroup <- colData %>% 
  mutate(FBCtrl_state_bin = ifelse(grepl('FB_Ctrl', condition), 1, 0)) %>% 
  select(5)


# binarize categorical variables

colData$tissue <- factor(colData$tissue, levels = c("FB","FG","HG","MT","MG","SG"))
colData$group <- factor(colData$group, levels = c("Ctrl","G2r","GFaFeces","GFa"))
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

group.out <- binarizeCategoricalColumns(colData$group,
                                         includePairwise = FALSE,
                                         includeLevelVsAll = TRUE,
                                         minCount = 1)

condition.out <- binarizeCategoricalColumns(colData$condition,
                                        includePairwise = FALSE,
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)

traits <- cbind(traits, tissue.out)
traits_group <- cbind(traits_group, group.out)
traits_tisgroup <- cbind(traits_tisgroup, condition.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

module.group.corr <- cor(module_eigengenes, traits_group, use = 'p')
module.group.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

module.traits_tg.corr <- cor(module_eigengenes, traits_tisgroup, use = 'p')
module.traits_tg.corr.pvals <- corPvalueStudent(module.traits_tg.corr, nSamples)



# visualize module-trait association as a heatmap

# TISSUES
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[39:44],
             y = names(heatmap.data)[1:38],
             col = c("blue1", "skyblue", "white", "pink", "red"))


module.gene.mapping <- as.data.frame(bwnet$colors)


# Fat bodies positively and negatively correlated modules (abs > 60)
FBpos_ye <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()
FBpos_dm <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkmagenta') %>% 
  rownames()
FBpos_dg <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkgreen') %>% 
  rownames()
FBpos_sb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'skyblue') %>% 
  rownames()
FBpos_rb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'royalblue') %>% 
  rownames()
FBpos_or <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'orange') %>% 
  rownames()
FBpos_wh <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'white') %>% 
  rownames()
FBneg_pi <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'pink') %>% 
  rownames()

# Fore gut positively and negatively correlated modules (abs > 60)
FGpos_gr <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'green') %>% 
  rownames()
FGpos_ta <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'tan') %>% 
  rownames()
FGneg_s3 <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'skyblue3') %>% 
  rownames()

# Mid gut positively correlated modules (abs > 60)
MGpos_bl <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'blue') %>% 
  rownames()
MGpos_dg <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkolivegreen') %>% 
  rownames()

# Hind gut positively correlated modules (abs > 60)
HGpos_pu <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'purple') %>% 
  rownames()
HGpos_mb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'midnightblue') %>% 
  rownames()
HGpos_ly <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'lightyellow') %>% 
  rownames()

# Malpighian tubes positively correlated modules (abs > 60)
MTpos_s3 <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'skyblue3') %>% 
  rownames()
MTpos_br <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()

# Salivary glans positively correlated modules (abs > 60)
SGpos_re <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames()
SGpos_dr <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkred') %>% 
  rownames()
SGpos_dg <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkolivegreen') %>% 
  rownames()

Module_list <- list(FBpos_ye,FBpos_dm,FBpos_dg,FBpos_sb,FBpos_rb,FBpos_or,FBpos_wh,FBneg_pi,
                    FGpos_gr,FGpos_ta,FGneg_s3,
                    MGpos_bl,MGpos_dg,
                    HGpos_pu,HGpos_mb,HGpos_ly,
                    MTpos_s3,MTpos_br,
                    SGpos_re,SGpos_dr,SGpos_dg)

Module_name <- c("FBpos_ye","FBpos_dm","FBpos_dg","FBpos_sb","FBpos_rb","FBpos_or","FBpos_wh","FBneg_pi",
                 "FGpos_gr","FGpos_ta","FGneg_s3",
                 "MGpos_bl","MGpos_dg",
                 "HGpos_pu","HGpos_mb","HGpos_ly",
                 "MTpos_s3","MTpos_br",
                 "SGpos_re","SGpos_dr","SGpos_dg")

# Codigo de 18.12.2025
for (i in 1:length(Module_list)) {
  modulo <- Module_list[[i]]
  modulo <- modulo[grepl("_", modulo)]
  if (length(modulo) >= 1) { 
    print(Module_name[[i]])
    print(length(Module_list[[i]]))
    print(paste(modulo, collapse = "|"))
  }
}

# [1] "FBpos_ye"
# [1] "CDS_defensin_g1_i1"                 "CDS_defensin_g1_i2"                 "CDS_defensin_g2"                   
# [4] "CDS_defensin_g8"                    "CDS_termicin_g2"                    "CDS_termicin_g3"                   
# [7] "CDS_blattellicin_g3"                "CDS_blattellicin_g4"                "IMD_SkpA_BgeEVm015521t1"           
# [10] "Blattellicin_g4_alt_BgeEVm012378t1" "Termicin_g1_alt_BgeEVm021689t1"     "Toll_Cact_BgeEVm007092t2"          
# [13] "Toll_spz6_BgeEVm007766t1"           "IMD_nemo_BgeEVm007246t1"            "IMD_PGRP-LC_BgeEVm007855t6"        
# [16] "Toll_Cact_BgeEVm007092t4"           "Toll_Cact_BgeEVm007092t5"           "Defensin_g2_alt_BgeEVm024201t2"    
# [19] "IMD_Duox_BgeEVm000625t3"            "IMD_Duox_BgeEVm000625t2"            "IMD_PGRP-LC_BgeEVm007855t3"        
# [22] "Toll_smt3_BgeEVm023921t1"  

# [1] "FBpos_dg"
# [1] "Toll_Cact_BgeEVm007092t3"

# [1] "FBpos_rb"
# [1] "Toll_Tl-Toll-like_BgeEVm001410t1" "Toll_Lwr_BgeEVm015656t2"          "Toll_Lwr_BgeEVm015656t1"         
# [4] "Toll_Tl-Toll-like_BgeEVm001410t2"

# [1] "FBpos_or"
# [1] "Toll_Ush_BgeEVm001356t1"

# [1] "FBneg_pi"
# [1] "CDS_defensin_g11"           "CDS_drosomycin_g5"          "IMD_akirin_BgeEVm014785t1" 
# [4] "Toll_cactin_BgeEVm003343t1" "Toll_tollip_BgeEVm011590t3" "Toll_tamo_BgeEVm003149t1"  
# [7] "IMD_Cyld_BgeEVm001283t3"   

# [1] "FGpos_gr"
# [1] "Toll_pli2_BgeEVm007692t1"         "Toll_Toll-like-R7_BgeEVm000745t1" "Toll_Tollo_BgeEVm000724t1"       
# [4] "Toll_Toll-like-R7_BgeEVm000745t3" "IMD_PGRP-LC_BgeEVm007855t10"      "IMD_PGRP-LC_BgeEVm007855t9"      
# [7] "Toll_Toll-like-R7_BgeEVm000745t2"

# [1] "FGpos_ta"
# [1] "Toll_dl_BgeEVm003323t2"     "Toll_Toll-6_BgeEVm000951t2" "IMD_PGRP-LC_BgeEVm007855t2"

# [1] "FGneg_s3"
# [1] "Toll_Tollo_BgeEVm000724t2"

# [1] "MGpos_bl"
# [1] "CDS_defensin_g9"                   "CDS_defensin_g10"                  "IMD_Skp2_BgeEVm010481t1"          
# [4] "IMD_eff_BgeEVm016400t1"            "Toll_Toll-like-R13_BgeEVm001898t3" "IMD_CASP_BgeEVm010772t1"          
# [7] "IMD_FADD2_BgeEVm012585t1"          "IMD_Npc2_BgeEVm016309t1"           "Toll_cactin_BgeEVm003343t2"       
# [10] "Toll_Toll-like-R13_BgeEVm001898t2" "IMD_CASP_BgeEVm010772t2"           "Toll_dl_BgeEVm003323t3"           
# [13] "Toll_MyD88_BgeEVm008366t3"         "IMD_Npc2_BgeEVm016309t2"           "Defensin_g10_alt_BgeEVm025224t2"  
# [16] "Drosomycin_g5_alt_BgeEVm035060t2" 

# [1] "HGpos_pu"
# [1] "Toll_Toll-6_BgeEVm000951t1" "Toll_spz3_BgeEVm007398t2"   "IMD_PGRP-LC_BgeEVm007855t5"

# [1] "HGpos_mb"
# [1] "CDS_drosomycin_g9"                "CDS_Attacin-like_g3A"             "Drosomycin_g9_alt_BgeEVm035060t1"
# [4] "IMD_rel_BgeEVm001701t2"      

# [1] "MTpos_s3"
# [1] "Toll_Tollo_BgeEVm000724t2"

# [1] "MTpos_br"
# [1] "CDS_drosomycin_g13"                "IMD_cad_BgeEVm010984t2"            "Toll_spz3_BgeEVm007398t1"         
# [4] "IMD_Ntf2_BgeEVm017640t1"           "IMD_PGRP-LC_BgeEVm007855t8"        "IMD_Cyld_BgeEVm001283t2"          
# [7] "IMD_cad_BgeEVm010984t1"            "IMD_Tab2_BgeEVm007139t3"           "IMD_Tab2_BgeEVm007139t2"          
# [10] "IMD_scny_BgeEVm001097t2"           "IMD_Ntf2_BgeEVm017640t2"           "Toll_tollip_BgeEVm011590t2"       
# [13] "Drosomycin_g13_alt_BgeEVm040569t1"

# [1] "SGpos_re"
# [1] "CDS_defensin_g3"                  "CDS_defensin_g4"                  "CDS_defensin_g5"                 
# [4] "CDS_defensin_g6"                  "CDS_defensin_g14"                 "CDS_defensin_g15"                
# [7] "CDS_defensin_g16_i1"              "CDS_defensin_g16_i2"              "IMD_PGRP-LC_BgeEVm007855t1"      
# [10] "Toll_Ush_BgeEVm001356t2"          "Toll_Aos1_BgeEVm009776t2"         "Defensin_g3-5_alt_BgeEVm027919t1"
# [13] "Defensin_g15_alt_BgeEVm023384t1" 

# GROUPS
heatmap.data <- merge(module_eigengenes, traits_group, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[39:42],
             y = names(heatmap.data)[1:38],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# Positively correlated with GFaFeces infected (all samples)
G2r_pos_pt <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'paleturquoise') %>% 
  rownames()
G2r_pos_sb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'skyblue') %>% 
  rownames()

Module_list <- list(G2r_pos_pt, G2r_pos_sb)
Module_name <- list("G2r_pos_pt", "G2r_pos_sb")

# Codigo de 18.12.2025
for (i in 1:length(Module_list)) {
  modulo <- Module_list[[i]]
  modulo <- modulo[grepl("_", modulo)]
  if (length(modulo) >= 1) { 
    print(Module_name[[i]])
    print(length(Module_list[[i]]))
    print(paste(modulo, collapse = "|"))
  }
}
# [1] "G2r_pos_pt"
# [1] 111
# [1] "IMD_Duox_BgeEVm000625t4"


# TISSUES + GROUP
heatmap.data.tg <- merge(module_eigengenes, traits_tisgroup, by = 'row.names')

head(heatmap.data.tg)

heatmap.data.tg <- heatmap.data.tg %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data.tg,
             x = names(heatmap.data.tg)[39:56],
             y = names(heatmap.data.tg)[1:38],
             col = c("blue1", "skyblue", "white", "pink", "red"))


# Positively and negatively correlated (all samples)
FBCtrl_pos_ma <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'darkmagenta') %>% 
  rownames()
FBCtrl_neg_vi <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'violet') %>% 
  rownames()
HGCtrl_pos_sa <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'salmon') %>% 
  rownames()
HGCtrl_pos_sb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'saddlebrown') %>% 
  rownames()
MTCtrl_pos_tu <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()
FBG2r_pos_pt <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'paleturquoise') %>% 
  rownames()
FBG2r_pos_sb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'skyblue') %>% 
  rownames()
HGFeces_pos_mb <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'midnightblue') %>% 
  rownames()
FBGF_pos_ma <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'magenta') %>% 
  rownames()

Module_list <- list(FBCtrl_pos_ma, FBCtrl_neg_vi, 
                    HGCtrl_pos_sa, HGCtrl_pos_sb,
                    MTCtrl_pos_tu, 
                    FBG2r_pos_pt, FBG2r_pos_sb,
                    HGFeces_pos_mb, 
                    FBGF_pos_ma)

Module_name <- list("FBCtrl_pos_ma", "FBCtrl_neg_vi",
                    "HGCtrl_pos_sa","HGCtrl_pos_sb",
                    "MTCtrl_pos_tu",
                    "FBG2r_pos_pt","FBG2r_pos_sb",
                    "HGFeces_pos_mb",
                    "FBGF_pos_ma")


# Codigo de 18.12.2025
for (i in 1:length(Module_list)) {
  modulo <- Module_list[[i]]
  modulo <- modulo[grepl("_", modulo)]
  if (length(modulo) >= 1) { 
    print(Module_name[[i]])
    print(length(Module_list[[i]]))
    print(paste(modulo, collapse = "|"))
  }
}

# [1] "MTCtrl_pos_tu"
# [1] 6700
# [1] "Toll_Lwr2_BgeEVm015620t1"
# [1] "FBG2r_pos_pt"
# [1] 111
# [1] "IMD_Duox_BgeEVm000625t4"
# [1] "HGFeces_pos_mb"
# [1] 700
# [1] "CDS_drosomycin_g9|CDS_Attacin-like_g3A|Drosomycin_g9_alt_BgeEVm035060t1|IMD_rel_BgeEVm001701t2"
# [1] "FBGF_pos_ma"
# [1] 1974
# [1] "CDS_defensin_g7|CDS_drosomycin_g3|CDS_Attacin-like_g2|CDS_Attacin-like_g3B|Defensin_g1_alt_BgeEVm024201t1|IMD_dsp1_BgeEVm009862t1|IMD_scny_BgeEVm001097t1|Toll_Ulp1_BgeEVm003986t1|Termicin_g3_alt_BgeEVm021689t2|Drosomcin_g3_alt_BgeEVm027887t1"

# Adding relevant Eggnog information of each module
Annotations <- read.csv2("../7.EggNOG/Bge_Annotations.tsv", sep = "\t")

Annotations_FBneg_dg <- Annotations[Annotations$query %in% FBneg_dg, ]
Annotations_FBneg_dg <- Annotations_FBneg_dg[,c("query", "Description", "Preferred_name", "PFAMs")]


# 6B. Intramodular analysis: Identifying driver genes ---------------


# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


MMcors <- as.data.frame(module.membership.measure)

# Calculate the gene significance and associated p-values

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
eggNOG <- read_tsv("../7.EggNOG/Bge_Annotations.tsv") %>%
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
ontology <- get_ontology(file = "../7.EggNOG/go.obo",
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
write_tsv(x = eggNOG, file = "../7.EggNOG/term2gene_GO.tsv")
write_tsv(x = eggNOG_term, file = "../7.EggNOG/term2name_GO.tsv")


# perform ORA
term2gene <- read_tsv("../7.EggNOG/term2gene_GO.tsv")
term2name <- read_tsv("../7.EggNOG/term2name_GO.tsv")


# KEGG ORA
eggNOG_kegg <- read_tsv("../7.EggNOG/Bge_Annotations.tsv") %>%
  dplyr::select(KEGG_ko, query) %>%
  dplyr::filter(KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  dplyr::mutate(gene = gsub("\\..*", "", query)) %>%
  dplyr::mutate(term = gsub("ko:", "", KEGG_ko)) %>%
  dplyr::select(term, gene) %>%
  distinct() %>%
  drop_na()




for(i in 1:length(Module_list)){ 
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  only_set = Module_list[[i]]
  title = Module_name[[i]]

  ############################################################################
  # ORA
  # read the gene list of interest
  over_enrichment <- enricher(only_set,
                              TERM2GENE = term2gene,
                              TERM2NAME = term2name,
                              pvalueCutoff = 0.05,
                              # universe = background_genes,
                              qvalueCutoff = 0.05)
  
  #save the enrichment result
  write.csv(file = paste0(title,"_ORA.csv"), x = over_enrichment@result)
  
  if (any(over_enrichment@result$p.adjust <= 0.05)){
    p <- dotplot(over_enrichment,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
      ggtitle("dotplot for GO ORA")
    
    ggsave(filename = paste0(title,"_ORA_dotplot.pdf"),            
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  }
  
  # Filtering
  over_enrichment_immuno <- over_enrichment@result
  over_enrichment_immuno <- over_enrichment_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd", 
                                                         over_enrichment_immuno$Description, 
                                                         ignore.case = TRUE),]
  over_enrichment_immuno <- over_enrichment_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell", 
                                                          over_enrichment_immuno$Description, 
                                                          ignore.case = TRUE),]
  over_enrichment_immuno <- over_enrichment_immuno %>% filter(p.adjust < 0.05)
  
  if (nrow(over_enrichment_immuno) > 0) {
    #save the enrichment result
    write.csv(file = paste0(title,"_Immuno_ORA.csv"), x = over_enrichment_immuno)
  }
  
  ############################################################################
  # KEGG
  # read the gene list of interest
  # create a list of kegg ortholog that I am interested in
  # over_set_kegg <- eggNOG_kegg %>%
  #   dplyr::filter(gene %in% only_set) %>%
  #   unlist() %>%
  #   as.vector()
  # 
  # # create a list of kegg ortholog that includes all kegg orthologs which form my background
  # background_kegg <- eggNOG_kegg$term
  # 
  # over_enrichment_kegg <- enrichKEGG(over_set_kegg,
  #                                    organism = "ko",
  #                                    keyType = "kegg",
  #                                    pvalueCutoff = 0.05,
  #                                    pAdjustMethod = "BH",
  #                                    universe = background_kegg,
  #                                    minGSSize = 10,
  #                                    maxGSSize = 500,
  #                                    qvalueCutoff = 0.05,
  #                                    use_internal_data = FALSE)
  # #save the enrichment result
  # write.csv(file = paste0(title,"_enrichment_KEGG.csv"),                
  #           x = over_enrichment_kegg@result)
  # 
  # if (any(over_enrichment_kegg@result$p.adjust <= 0.05)){
  #   p <- dotplot(over_enrichment_kegg,
  #                x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
  #                color="p.adjust",
  #                orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
  #                showCategory=100,
  #                font.size=8) +
  #     ggtitle("dotplot for KEGG ORA")
  #   
  #   ggsave(filename = paste0(title,"_enrichment_KEGG.pdf"),                
  #          plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
  # }
}



