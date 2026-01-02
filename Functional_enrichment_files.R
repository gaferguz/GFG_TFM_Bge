#########################
# FUNCTIONAL ENRICHMENT #
#########################

## Enriquecimiento funcional usando eggnog

install.packages(c("BiocManager", "devtools"))
library("devtools")
library("BiocManager")
BiocManager::install("tidyverse")
devtools::install_github("GuangchuangYu/GOSemSim")
devtools::install_github("GuangchuangYu/clusterProfiler")

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



