#scripts for deming regression model building
library(tidyverse) 
library(deming)
library(ggplot2)
library(tidyverse)
library(igraph)
library(metafor)

#create full_gse_list
#top_genes_WT2_gse102348 - outlier
names_vec = c()
full_gse_list = list()
for (i in list.files("full_FClists/")){
  print(i)
  flops = readRDS(paste("full_FClists/", i, sep=""))
  full_gse_list = append(full_gse_list, list(flops))
  names_vec = c(names_vec, i)
}
names(full_gse_list) = names_vec

for (i in full_gse_list){
  rownames_vec = union(rownames(i), rownames_vec)
}
colnames_vec = names(full_gse_list)

logFCmatrixregr = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(logFCmatrixregr) = rownames_vec
colnames(logFCmatrixregr) = colnames_vec

SEmatrixregr = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(SEmatrixregr) = rownames_vec
colnames(SEmatrixregr) = colnames_vec

for (name in names(full_gse_list)){
  logFCmatrixregr[rownames(full_gse_list[[name]]), name] =  full_gse_list[[name]]$logFC
  SEmatrixregr[rownames(full_gse_list[[name]]), name] =  full_gse_list[[name]]$SE
}

# normalize by sd:
for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}

# make totalrownamematrix:
totalrownamematrix = list()
for (i in 1:(length(full_gse_list)-1)){
  for (j in (i + 1):length(full_gse_list)){
    topA = full_gse_list[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = full_gse_list[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    rownamesA = rownames(subset(full_gse_list[[i]], rownames(full_gse_list[[i]]) %in% totalrownames))
    rownamesB = rownames(subset(full_gse_list[[j]], rownames(full_gse_list[[j]]) %in% totalrownames))
    totalrownamematrix[[names(full_gse_list)[i]]][[names(full_gse_list)[j]]] = intersect(rownamesA, rownamesB)
    totalrownamematrix[[names(full_gse_list)[j]]][[names(full_gse_list)[i]]] = intersect(rownamesA, rownamesB)
  }
}

cormethod = "pearson"
corpval750 = data.frame()
cormatrix750 = data.frame()
coradjpval750 = data.frame()
thres = "750"

#CORMATRIX
for (i in 1:length(full_gse_list)){
  for (j in i:length(full_gse_list)){
    topA = full_gse_list[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = full_gse_list[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    cormatrix750[names(full_gse_list)[i], names(full_gse_list)[j]] = cor(full_gse_list[[i]][totalrownames,]$logFC, full_gse_list[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    cormatrix750[names(full_gse_list)[j], names(full_gse_list)[i]] = cor(full_gse_list[[i]][totalrownames,]$logFC, full_gse_list[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
  }
}


#CORPVAL
for (i in 1:length(full_gse_list)){
  for (j in i:length(full_gse_list)){
    topA = full_gse_list[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = full_gse_list[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    corpval750[names(full_gse_list)[i], names(full_gse_list)[j]] = as.numeric(cor.test(full_gse_list[[i]][totalrownames,]$logFC, full_gse_list[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    corpval750[names(full_gse_list)[j], names(full_gse_list)[i]] = as.numeric(cor.test(full_gse_list[[i]][totalrownames,]$logFC, full_gse_list[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
  }
}

#CORADJPVAL
vec = as.vector(corpval750[upper.tri(corpval750, diag = F)])
vec = p.adjust(vec, method = "BH")
tempmatrix = matrix(0, length(full_gse_list), length(full_gse_list))
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

saveRDS(logFCmatrixregr, file = "logFCmatrixregr_full.Rds")
saveRDS(SEmatrixregr, file = "SEmatrixregr_full.Rds")
saveRDS(cortestsign750, file = "cortestsign750_full.Rds")
saveRDS(totalrownamematrix, file = "totalrownamematrix_full.Rds")

# plot exapmles:
kekmatrix = cormatrix750 %>% rownames_to_column("datasetid1")
kekmatrix = gather(kekmatrix, datasetid2, corvalue, -datasetid1)
kekmatrix$corvalue = as.numeric(as.character(kekmatrix$corvalue))
kekmatrix = kekmatrix %>% filter(corvalue != 1) %>% top_n(10, corvalue) %>% distinct(corvalue, .keep_all = T)

for (i in 1:length(rownames(kekmatrix))){
  plot1 = deming(logFCmatrix_full_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid2"]] ~ logFCmatrix_full_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid1"]] - 1)
  
  ggheatmap = ggplot(logFCmatrix_full_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],], aes_string(x = kekmatrix[i,"datasetid1"], y = kekmatrix[i,"datasetid2"])) + geom_point() +
    geom_abline(slope = plot1$coefficients[2], intercept = 0, colour = "blue", size = 1) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_abline(slope = kres[which(colnames(logFCmatrix_full_goodboys) == kekmatrix[i, "datasetid2"])]/kres[which(colnames(logFCmatrix_full_goodboys) == kekmatrix[i, "datasetid1"])], intercept = 0, colour = "red", size = 1)
  print(ggheatmap)
  pdf(paste0("demingexample_all_", i,  ".pdf"))
  print(ggheatmap)
  dev.off()
}

deming_result_full = readRDS("deming_result_full.Rds")
logFCmatrixregr_full = readRDS("logFCmatrixregr_full.Rds")
SEmatrixregr_full = readRDS("SEmatrixregr_full.Rds")

kres = deming_result_full
for (i in 1:length(colnames(logFCmatrixregr_full))){
  SEmatrixregr_full[,i] = SEmatrixregr_full[,i] / kres[i]
  logFCmatrixregr_full[,i] = logFCmatrixregr_full[,i] / kres[i]
}
temp_colnames = colnames(logFCmatrixregr_full)
logFCmatrixregr_full = data.frame(logFCmatrixregr_full)
logFCmatrixregr_full$NACount = rowSums(is.na(logFCmatrixregr_full))
ggplot(logFCmatrixregr_full, aes(x = NACount)) + geom_density()
logFCmatrix_full_goodboys = subset(logFCmatrixregr_full, logFCmatrixregr_full$NACount <=round((2*ncol(logFCmatrixregr_full))/3))
SEmatrix_full_goodboys = subset(SEmatrixregr_full, logFCmatrixregr_full$NACount <=round((2*ncol(logFCmatrixregr_full))/3))
logFCmatrix_full_goodboys$NACount = NULL
colnames(logFCmatrix_full_goodboys) = temp_colnames

sourcedata_full = as.data.frame(colnames(logFCmatrix_full_goodboys))
rownames(sourcedata_full) = colnames(logFCmatrix_full_goodboys)
colnames(sourcedata_full) = "datasets"
cops = c()
for (i in strsplit(rownames(sourcedata_full), "_")){
  cops = c(cops, gsub(".Rds", "", i[length(i)]))
}
sourcedata_full[,"datasets"] = cops

full_signature = signature_builder(logFCmatrix_full_goodboys, SEmatrix_full_goodboys, sourcedata_full)
saveRDS(full_signature, file = "full_signature.Rds")
