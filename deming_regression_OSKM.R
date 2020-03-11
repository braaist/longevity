#scripts for deming regression model building
library(tidyverse) 
library(deming)
library(ggplot2)
library(tidyverse)
library(igraph)
library(metafor)

#create OSKM_gse_list
#top_genes_HEP_114581 - outlier
names_vec = c()
OSKM_gse_list = list()
for (i in list.files("OSKM_FClists/")){
  print(i)
  flops = readRDS(paste("OSKM_FClists/", i, sep=""))
  OSKM_gse_list = append(OSKM_gse_list, list(flops))
  names_vec = c(names_vec, i)
}
names(OSKM_gse_list) = names_vec

for (i in OSKM_gse_list){
  rownames_vec = union(rownames(i), rownames_vec)
}
colnames_vec = names(OSKM_gse_list)

logFCmatrixregr = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(logFCmatrixregr) = rownames_vec
colnames(logFCmatrixregr) = colnames_vec

SEmatrixregr = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(SEmatrixregr) = rownames_vec
colnames(SEmatrixregr) = colnames_vec

for (name in names(OSKM_gse_list)){
  logFCmatrixregr[rownames(OSKM_gse_list[[name]]), name] =  OSKM_gse_list[[name]]$logFC
  SEmatrixregr[rownames(OSKM_gse_list[[name]]), name] =  OSKM_gse_list[[name]]$SE
}

# normalize by sd:
for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}

# make totalrownamematrix:
totalrownamematrix = list()
for (i in 1:(length(OSKM_gse_list)-1)){
  for (j in (i + 1):length(OSKM_gse_list)){
    topA = OSKM_gse_list[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = OSKM_gse_list[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    rownamesA = rownames(subset(OSKM_gse_list[[i]], rownames(OSKM_gse_list[[i]]) %in% totalrownames))
    rownamesB = rownames(subset(OSKM_gse_list[[j]], rownames(OSKM_gse_list[[j]]) %in% totalrownames))
    totalrownamematrix[[names(OSKM_gse_list)[i]]][[names(OSKM_gse_list)[j]]] = intersect(rownamesA, rownamesB)
    totalrownamematrix[[names(OSKM_gse_list)[j]]][[names(OSKM_gse_list)[i]]] = intersect(rownamesA, rownamesB)
  }
}

cormethod = "pearson"
corpval750 = data.frame()
cormatrix750 = data.frame()
coradjpval750 = data.frame()
thres = "750"

#CORMATRIX
for (i in 1:length(OSKM_gse_list)){
  for (j in i:length(OSKM_gse_list)){
    topA = OSKM_gse_list[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = OSKM_gse_list[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    cormatrix750[names(OSKM_gse_list)[i], names(OSKM_gse_list)[j]] = cor(OSKM_gse_list[[i]][totalrownames,]$logFC, OSKM_gse_list[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    cormatrix750[names(OSKM_gse_list)[j], names(OSKM_gse_list)[i]] = cor(OSKM_gse_list[[i]][totalrownames,]$logFC, OSKM_gse_list[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
  }
}


#CORPVAL
for (i in 1:length(OSKM_gse_list)){
  for (j in i:length(OSKM_gse_list)){
    topA = OSKM_gse_list[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = OSKM_gse_list[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    corpval750[names(OSKM_gse_list)[i], names(OSKM_gse_list)[j]] = as.numeric(cor.test(OSKM_gse_list[[i]][totalrownames,]$logFC, OSKM_gse_list[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    corpval750[names(OSKM_gse_list)[j], names(OSKM_gse_list)[i]] = as.numeric(cor.test(OSKM_gse_list[[i]][totalrownames,]$logFC, OSKM_gse_list[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
  }
}

#CORADJPVAL
vec = as.vector(corpval750[upper.tri(corpval750, diag = F)])
vec = p.adjust(vec, method = "BH")
tempmatrix = matrix(0, length(OSKM_gse_list), length(OSKM_gse_list))
tempmatrix[lower.tri(tempmatrix, diag = F)] = vec
tempmatrix = t(tempmatrix)
tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
coradjpval750 = tempmatrix
colnames(coradjpval750) = colnames(corpval750)
rownames(coradjpval750) = rownames(corpval750)


cortestsign750 = data.frame()
for (colname in colnames(cormatrix750)){
  for (rowname in rownames(cormatrix750)){
    if (coradjpval750[rowname, colname] < 0.05){
      if (cormatrix750[rowname, colname] > 0.1){
        cortestsign750[rowname, colname] = 1
      } else if (cormatrix750[rowname, colname] < -0.1){
        cortestsign750[rowname, colname] = -1
      } else {
        cortestsign750[rowname, colname] = 0
      }
    } else{
      cortestsign750[rowname, colname] = 0
    }
  }
}


#make graph
matrixfornetwork = cortestsign750
for (i in 1:length(rownames(matrixfornetwork))){
  for (j in 1:length(colnames(matrixfornetwork))){
    if (matrixfornetwork[i, j] == -1){
      matrixfornetwork[i, j] = 0
    }
  }
}
network = graph_from_adjacency_matrix(as.matrix(abs(matrixfornetwork)), mode = "undirected", diag = F)
plot(network, vertex.size = 2, vertex.label.cex = 0.6, edge.width = 1)

saveRDS(logFCmatrixregr, file = "logFCmatrixregr_OSKM.Rds")
saveRDS(SEmatrixregr, file = "SEmatrixregr_OSKM.Rds")
saveRDS(cortestsign750, file = "cortestsign750_OSKM.Rds")
saveRDS(totalrownamematrix, file = "totalrownamematrix_OSKM.Rds")

# plot exapmles:
kekmatrix = cormatrix750 %>% rownames_to_column("datasetid1")
kekmatrix = gather(kekmatrix, datasetid2, corvalue, -datasetid1)
kekmatrix$corvalue = as.numeric(as.character(kekmatrix$corvalue))
kekmatrix = kekmatrix %>% filter(corvalue != 1) %>% top_n(10, corvalue) %>% distinct(corvalue, .keep_all = T)

for (i in 1:length(rownames(kekmatrix))){
  plot1 = deming(logFCmatrix_OSKM_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid2"]] ~ logFCmatrix_OSKM_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid1"]] - 1)
  
  ggheatmap = ggplot(logFCmatrix_OSKM_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],], aes_string(x = kekmatrix[i,"datasetid1"], y = kekmatrix[i,"datasetid2"])) + geom_point() +
    geom_abline(slope = plot1$coefficients[2], intercept = 0, colour = "blue", size = 1) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_abline(slope = kres[which(colnames(logFCmatrix_OSKM_goodboys) == kekmatrix[i, "datasetid2"])]/kres[which(colnames(logFCmatrix_OSKM_goodboys) == kekmatrix[i, "datasetid1"])], intercept = 0, colour = "red", size = 1)
  print(ggheatmap)
  pdf(paste0("demingexample_OSKM_", i,  ".pdf"))
  print(ggheatmap)
  dev.off()
}


deming_result_OSKM = readRDS("deming_result_OSKM.Rds")
logFCmatrixregr_OSKM = readRDS("logFCmatrixregr_OSKM.Rds")
SEmatrixregr_OSKM = readRDS("SEmatrixregr_OSKM.Rds")

kres = deming_result_OSKM
for (i in 1:length(colnames(logFCmatrixregr_OSKM))){
  SEmatrixregr_OSKM[,i] = SEmatrixregr_OSKM[,i] / kres[i]
  logFCmatrixregr_OSKM[,i] = logFCmatrixregr_OSKM[,i] / kres[i]
}
logFCmatrixregr_OSKM = data.frame(logFCmatrixregr_OSKM)
logFCmatrixregr_OSKM = data.frame(logFCmatrixregr_OSKM)
logFCmatrixregr_OSKM$NACount = rowSums(is.na(logFCmatrixregr_OSKM))
ggplot(logFCmatrixregr_OSKM, aes(x = NACount)) + geom_density()
logFCmatrix_OSKM_goodboys = subset(logFCmatrixregr_OSKM, logFCmatrixregr_OSKM$NACount <=round((2*ncol(logFCmatrixregr_OSKM))/3))
SEmatrix_OSKM_goodboys = subset(SEmatrixregr_OSKM, logFCmatrixregr_OSKM$NACount <=round((2*ncol(logFCmatrixregr_OSKM))/3))
logFCmatrix_OSKM_goodboys$NACount = NULL
colnames(logFCmatrix_OSKM_goodboys) = temp_colnames

sourcedata_OSKM = as.data.frame(colnames(logFCmatrix_OSKM_goodboys))
rownames(sourcedata_OSKM) = colnames(logFCmatrix_OSKM_goodboys)
colnames(sourcedata_OSKM) = "datasets"
cops = c()
for (i in strsplit(rownames(sourcedata_OSKM), "_")){
  cops = c(cops, gsub(".Rds", "", i[length(i)]))
}
sourcedata_OSKM[,"datasets"] = cops
OSKM_signature = signature_builder(logFCmatrix_OSKM_goodboys, SEmatrix_OSKM_goodboys, sourcedata_OSKM)
saveRDS(OSKM_signature, file = "OSKM_signature.Rds")



