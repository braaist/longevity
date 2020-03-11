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


#SASHA's scripts
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

#SASHA's functions for cormat
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#Stan's scripts
#Deming regression minimizer
deming_minimizer = function(logFCmatrixregr){
  fn = function(k_no_first){
    k = c()
    k[1] = 1
    k[2:length(colnames(logFCmatrixregr))] = k_no_first
    res = 0
    for (i in 1:(length(colnames(logFCmatrixregr))-1)){
      namei = colnames(logFCmatrixregr)[i]
      for (j in (i + 1):length(colnames(logFCmatrixregr))){
        namej = colnames(logFCmatrixregr)[j]
        if (cortestsign500[namei, namej] != 1){
          next
        }
        totalrownames = totalrownamematrix[[namei]][[namej]]
        ai = logFCmatrixregr[totalrownames, namei]
        aj = logFCmatrixregr[totalrownames, namej]
        res = res + sum(
          (((aj - (k[j]/k[i])*ai)^2)*((ai - (k[i]/k[j])*aj)^2))/
            (((aj - (k[j]/k[i])*ai)^2)+((ai - (k[i]/k[j])*aj)^2)))/length(totalrownames)
      }
    }
    return(res)
  }
  
  kvec = rnorm(length(colnames(logFCmatrixregr)) - 1, 1, 1)
  ptm <- proc.time()
  optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B")
  proc.time() - ptm
  kres = c(1, optimized$par)
  minimum = optimized$value
  bigres = list(kres, minimum)
  names(bigres) = c("coefs", "minimum")
  return(bigres)
}

# mixed effect model signature builder

signature_builder = function(logFCmatrixregr_sb, SEmatrixregr_sb, sourcedata_sb){
  goodgenes = c()
  signature = data.frame()
  genenumber = 0
  for (genename in rownames(logFCmatrixregr_sb)){
    genenumber = genenumber + 1
    percentready = (genenumber/length(rownames(logFCmatrixregr_sb))) * 100
    if (genenumber %% 1000 == 0){
      print(paste0("I'm on gene No. ", genenumber, " (", round(percentready, 2), "% done)"))
    }
    logFC = logFCmatrixregr_sb[genename,]
    logFC = logFC[!is.na(logFC)]
    SE = SEmatrixregr_sb[genename,]
    SE = SE[!is.na(SE)]
    sourcevec = as.factor(sourcedata_sb[colnames(logFCmatrixregr_sb)[!is.na(logFCmatrixregr_sb[genename,])],])

    tryCatch(
      {
        mixedeffres = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec))
        signature = rbind(signature, c(mixedeffres$b[1], mixedeffres$pval))
        goodgenes = c(goodgenes, genename)
      },
      error=function(cond) {
        message("Fucked up")
        message("Here's the original error message:")
        message(cond)
      }
    )
  }
  #rownames(signature) = totalgenes[-which(totalgenes %in% badgenes)]
  rownames(signature) = goodgenes
  colnames(signature) = c("logFC", "pval")
  
  signature$adj_pval = p.adjust(signature$pval, method = "BH")
  return(signature)
}


cormat_maker <- function(list_for_cormat){
  cormat = data.frame()
  for (i in names(list_for_cormat)){
    print(i)
    for (j in names(list_for_cormat)){
      pops = (union(rownames(list_for_cormat[[i]])[1:500], rownames(list_for_cormat[[j]])[1:500]))
      temp_m = data.frame()
      for (m in pops){
        temp_m[m,1] = list_for_cormat[[i]][m,1]
        temp_m[m,2] = list_for_cormat[[j]][m,1]
      }
      print(cor(x = temp_m[,1], y = temp_m[,2], use = "complete.obs", method = "spearman"))
      cormat[i,j] = cor(x = temp_m[,1], y = temp_m[,2], use = "complete.obs", method = "spearman")
    }
  }
  cormatrix = apply(cormat, 2, rev)
  upper_tri <- get_upper_tri(cormatrix)
  # Melt the correlation matrix
  #melted_cormat <- melt(upper_tri, na.rm = TRUE)
  melted_cormat <- melt(cormatrix, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  ggheatmap
}


cormat_maker_full <- function(list_for_cormat){
  cormat = data.frame()
  for (i in names(list_for_cormat)){
    print(i)
    for (j in names(list_for_cormat)){
      pops = (union(rownames(list_for_cormat[[i]]), rownames(list_for_cormat[[j]])))
      temp_m = data.frame()
      for (m in pops){
        temp_m[m,1] = list_for_cormat[[i]][m,1]
        temp_m[m,2] = list_for_cormat[[j]][m,1]
      }
      print(cor(x = temp_m[,1], y = temp_m[,2], use = "complete.obs", method = "spearman"))
      cormat[i,j] = cor(x = temp_m[,1], y = temp_m[,2], use = "complete.obs", method = "spearman")
    }
  }
  cormatrix = apply(cormat, 2, rev)
  upper_tri <- get_upper_tri(cormatrix)
  # Melt the correlation matrix
  #melted_cormat <- melt(upper_tri, na.rm = TRUE)
  melted_cormat <- melt(cormatrix, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  ggheatmap
}



signature_builder_2terms = function(logFCmatrixregr_sb, SEmatrixregr_sb, sourcedata_sb){
  goodgenes = c()
  signature = data.frame()
  genenumber = 0
  for (genename in rownames(logFCmatrixregr_sb)){
    genenumber = genenumber + 1
    percentready = (genenumber/length(rownames(logFCmatrixregr_sb))) * 100
    if (genenumber %% 1000 == 0){
      print(paste0("I'm on gene No. ", genenumber, " (", round(percentready, 2), "% done)"))
    }
    logFC = logFCmatrixregr_sb[genename,]
    logFC = logFC[!is.na(logFC)]
    SE = SEmatrixregr_sb[genename,]
    SE = SE[!is.na(SE)]
    sourcevec = as.factor(sourcedata_sb[colnames(logFCmatrixregr_sb)[!is.na(logFCmatrixregr_sb[genename,])],])
    
    tryCatch(
      {
        mixedeffres = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec))
        signature = rbind(signature, c(mixedeffres$b[1], mixedeffres$pval))
        goodgenes = c(goodgenes, genename)
      },
      error=function(cond) {
        message("Fucked up")
        message("Here's the original error message:")
        message(cond)
      }
    )
  }
  #rownames(signature) = totalgenes[-which(totalgenes %in% badgenes)]
  rownames(signature) = goodgenes
  colnames(signature) = c("logFC", "pval")
  
  signature$adj_pval = p.adjust(signature$pval, method = "BH")
  return(signature)
}




