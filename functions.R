#FUNCTIONS USED IN THIS PROJECT
#get GSE file create phenodata table
library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
library(biomaRt)
library(mogene20sttranscriptcluster.db)

entrez_converter <- function(ACCNUM){
  #function for convertation ACCNUM to ENTREZID
  if(is.null(xx[[as.character(ACCNUM)]])){
    return(NA)
  }else{
    return(xx[[as.character(ACCNUM)]])}
}

limma_maker <- function(assayData, des_matrix, cont_matrix){
  #function for getting data about differential expression and creating volcano plots
  fit = lmFit(assayData, des_matrix)
  fit_contrast = eBayes(contrasts.fit(fit, cont_matrix))
  volcanoplot(fit_contrast)
  top_genes = topTable(fit_contrast, number = Inf, adjust = "BH", confint = TRUE)
  top_genes$SE = (top_genes$CI.R-top_genes$CI.L)/3.92
  top_genes_limma <<- top_genes
  return(top_genes)
}

ensembl2entrez_mm <- function(ENSEMBL_IDS){
  ENTREZ_IDS = c()
  x <- org.Mm.egENSEMBL2EG
  mapped_probes = mappedkeys(x)
  xx = as.list(x[mapped_probes])
  for (i in ENSEMBL_IDS){
  if(is.null(xx[[as.character(i)]])){
    ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }
  else if(length(xx[[as.character(i)]]) != 1){
    ENTREZ_IDS = c(ENTREZ_IDS, NA)
  }else{
    ENTREZ_IDS = c(ENTREZ_IDS, xx[[as.character(i)]])
    }
  }
  return(ENTREZ_IDS)
}
affy2entrez_mm <- function(AFFY_IDS){
  ENTREZ_IDS = c()
  x <- mouse430a2ENTREZID
  mapped_probes = mappedkeys(x)
  xx = as.list(x[mapped_probes])
  for (i in AFFY_IDS){
    if(is.null(xx[[as.character(i)]])){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }
    else if(length(xx[[as.character(i)]]) != 1){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }else{
      ENTREZ_IDS = c(ENTREZ_IDS, xx[[as.character(i)]])
    }
  }
  return(ENTREZ_IDS)
}

ilmna2entrez_hs <- function(ILMN_IDS){
  ENTREZ_IDS = c()
  x <- illuminaHumanv4ENTREZID
  mapped_probes = mappedkeys(x)
  xx = as.list(x[mapped_probes])
  for (i in ILMN_IDS){
    if(is.null(xx[[as.character(i)]])){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }
    else if(length(xx[[as.character(i)]]) != 1){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }else{
      ENTREZ_IDS = c(ENTREZ_IDS, xx[[as.character(i)]])
    }
  }
  return(ENTREZ_IDS)
}

ilmnav32entrez_hs <- function(ILMN_IDS){
  ENTREZ_IDS = c()
  x <- illuminaHumanv3ENTREZID
  mapped_probes = mappedkeys(x)
  xx = as.list(x[mapped_probes])
  for (i in ILMN_IDS){
    if(is.null(xx[[as.character(i)]])){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }
    else if(length(xx[[as.character(i)]]) != 1){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }else{
      ENTREZ_IDS = c(ENTREZ_IDS, xx[[as.character(i)]])
    }
  }
  return(ENTREZ_IDS)
}

mogene2entrez_mm <- function(AFFY_IDS){
  ENTREZ_IDS = c()
  x <- mogene20sttranscriptclusterENTREZID
  mapped_probes = mappedkeys(x)
  xx = as.list(x[mapped_probes])
  for (i in AFFY_IDS){
    if(is.null(xx[[as.character(i)]])){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }
    else if(length(xx[[as.character(i)]]) != 1){
      ENTREZ_IDS = c(ENTREZ_IDS, NA)
    }else{
      ENTREZ_IDS = c(ENTREZ_IDS, xx[[as.character(i)]])
    }
  }
  return(ENTREZ_IDS)
}



