# Cargamos las librerias necesarias 
library(dplyr)
library(ggplot2)
library(data.table)
library(MatrixGenerics)
library(matrixStats)
library(ggpubr)


##################################
#           FUNCIONES            #
##################################

TPM_filtering <- function(filtroTPM, datos, factor, TPMs){
  # Inicalizacion de un vector que conteanga los ids de los transcritos
  Passed_tx <- c()
  # Iteracion usando como subgrupo la columna "factor" de "datos"
  for (group in unique(datos[,factor])) {
    subdatos <- datos[grepl(group, datos[,factor]),]
    submatrix <- TPMs[colnames(TPMs) %in% subdatos$sample]
    threshold <- round(3/4*(ncol(submatrix)),0)
    # Binarizacion de la matriz
    # Valores >= 1 seran 1, valores menores seran 0, luego hacemos rowsums
    # Iteramos columna a columna para estar seguros de la conversion numerica
    submatrix <- as.matrix(submatrix)
    class(submatrix) <- "numeric"
    submatrix[submatrix < filtroTPM] <- 0
    submatrix[submatrix >= filtroTPM] <- 1
    # Despues de binarizar sumamos las filas y si ésta es >= theshold, el transcrito pasa el filtro
    Sum <- rowSums(submatrix)
    subdf_sum <- data.frame(TPMs$target_id, Sum)
    subdf_sumf2 <- subdf_sum %>% filter(Sum >= threshold)
    colnames(subdf_sumf2) <- c("target_id","Sum")
    # Se almacenan los ids de los transcritos que pasan y se eliminan redundancias
    Passed_tx <- c(Passed_tx, subdf_sumf2$target_id)
    Passed_tx <- unique(Passed_tx)
    # Informacion que aparece por pantalla
    print(paste0("N de tx: ",length(Passed_tx), " - Grupo: ", group, "- Umbral 75: ", threshold))
  }
  return(Passed_tx)
}


########################
#         MAIN         #
########################

datos <- read.csv2("../samples.txt", sep = "\t")

factor <- c("condition")

# Esta misma funcion se uso con matrices de abundancia TPM sacadas del script de kallisto para:

# 1 : Primeros ensamblajes filtrados (para Trinity y Spades)

# 2 : Filtrado de los transcritos de EvidentialGene tr2aacds.pl (okay + cull)

# Ejemplo

TPMs <- read.csv2("../Kallisto/Kallisto_Results/Kallisto_Merged/Abundance.tsv", sep = " ")
# modificacion de las columnas para que sean iguales a los nombres de las muestras en los metadatos
colnames(TPMs) <- gsub("\\.","-",
                       gsub("tpm_","", 
                            gsub("_S.*","", colnames(TPMs))))


Transcritos_seleccionados <- TPM_filtering(1, datos, factor, TPMs)

# Escribimos un archivo de texto con los transcritos seleccionados para filtrarlos con seqtk
data.table::fwrite(list(Transcritos_seleccionados), file = "./Okay_over1TPM.txt")
