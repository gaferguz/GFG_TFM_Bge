# Leemos los transcritos cull y okay como dos tablas separadas
cull_meta <-  read.table("./cullset.tsv", header = FALSE, sep = "\t", )
colnames(cull_meta) <- c("Public_mRNA_ID","originalID","PublicGeneID","AltNum",	
                         "Class","AAqual","pIdAln","Notes","Oids")

# Para los cull solo nos interesaria las isoformas no descartadas
cull_meta <- cull_meta[!grepl("cullalt", cull_meta$Class),]

okay_meta <-  read.table("./okayset.tsv", header = FALSE, sep = "\t", )
colnames(okay_meta) <- c("Public_mRNA_ID","originalID","PublicGeneID","AltNum",	
                         "Class","AAqual","pIdAln","Notes","Oids")

# Extraemos los ids de los transcritos de interes
cull_tx <- cull_meta$Public_mRNA_ID
okay_tx <- okay_meta$Public_mRNA_ID

okaycullset <- c(okay_tx, cull_tx) 
length(unique(gsub("t.*","",over1TPM_okaycull))) # Ver los locus totales

data.table::fwrite(list(okaycullset), file = "./Okay_cullset.txt")
