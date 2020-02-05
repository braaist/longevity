#FUNCTIONS USED IN THIS PROJECT
#get GSE file create phenodata table
library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
library(biomaRt)
library(mogene20sttranscriptcluster.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


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

ensembl2entrez_hs <- function(ENSEMBL_IDS){
  ENTREZ_IDS = c()
  x <- org.Hs.egENSEMBL
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


#SASHA's script
#Download_raw_reads is a function for downloading RNA-seq reads from featureCounts output
#phenodata rownames should be the same as names of featureCounts count files

Download_raw_reads <- function(featurecounts_dir,phenodata){
  setwd(featurecounts_dir)
  
  temp_data <- read.csv(dir()[!grepl(".summary$",dir())][1],header=T,sep='\t',skip = 1)
  counts_star <- data.frame(ID=temp_data$Geneid)
  rownames(counts_star) <- counts_star$ID
  for (i in dir()[!grepl(".summary$",dir())]){
    temp_name <- strsplit(i,".count")[[1]][1]
    temp_name <- gsub("-","_",temp_name)
    temp_data <- read.csv(i,header=T,sep='\t',skip = 1)
    temp_counts <- temp_data[,7]
    counts_star[temp_name] <- temp_counts
  }
  rm(temp_data)
  counts_star$ID <- rownames(counts_star)
  counts_star <- counts_star[,-1]
  inter_samples <- intersect(rownames(phenodata),colnames(counts_star))
  counts_star <- counts_star[,inter_samples]
  
  #Create expression set object
  library("lumi")
  phenodata <- phenodata[inter_samples,]
  meta.info <- data.frame(labelDescription=colnames(phenodata))
  pheno <- new("AnnotatedDataFrame", data = phenodata, varMetadata = meta.info)
  RNAseq_counts_star <- new("ExpressionSet", exprs=as.matrix(counts_star),phenoData=pheno)
  
  return(RNAseq_counts_star)
}

