#scripts for deming regression model building in all DS
library(tidyverse) 
library(deming)
library(ggplot2)
library(tidyverse)
library(igraph)
library(metafor)


#first step

#TRESHOLD 750 WAS CHOSEN
#create gse_list_main
names_vec = c()
gse_list_main = list()
for (i in list.files("all_FClists2/")){
  print(i)
  flops = readRDS(paste("all_FClists2/", i, sep=""))
  gse_list_main = append(gse_list_main, list(flops))
  names_vec = c(names_vec, i)
}
names(gse_list_main) = names_vec

rownames_vec = c()
for (i in gse_list_main){
  rownames_vec = union(rownames(i), rownames_vec)
}
View(rownames_vec)
colnames_vec = names(gse_list_main)

SEmatrixregr = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(SEmatrixregr) = rownames_vec
colnames(SEmatrixregr) = colnames_vec

logFCmatrixregr = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(logFCmatrixregr) = rownames_vec
colnames(logFCmatrixregr) = colnames_vec

for (name in names( gse_list_main)){
  SEmatrixregr[rownames( gse_list_main[[name]]), name] =  gse_list_main[[name]]$SE
  logFCmatrixregr[rownames( gse_list_main[[name]]), name] =  gse_list_main[[name]]$logFC
}

# normalize by sd:
for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}

# make totalrownamematrix:
totalrownamematrix = list()
for (i in 1:(length(gse_list_main)-1)){
  for (j in (i + 1):length(gse_list_main)){
    topA = gse_list_main[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = gse_list_main[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    rownamesA = rownames(subset(gse_list_main[[i]], rownames(gse_list_main[[i]]) %in% totalrownames))
    rownamesB = rownames(subset(gse_list_main[[j]], rownames(gse_list_main[[j]]) %in% totalrownames))
    totalrownamematrix[[names(gse_list_main)[i]]][[names(gse_list_main)[j]]] = intersect(rownamesA, rownamesB)
    totalrownamematrix[[names(gse_list_main)[j]]][[names(gse_list_main)[i]]] = intersect(rownamesA, rownamesB)
  }
}

cormethod = "pearson"
corpval750 = data.frame()
cormatrix750 = data.frame()
coradjpval750 = data.frame()
thres = "750"

#CORMATRIX
for (i in 1:length(gse_list_main)){
  for (j in i:length(gse_list_main)){
    topA = gse_list_main[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = gse_list_main[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    cormatrix750[names(gse_list_main)[i], names(gse_list_main)[j]] = cor(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    cormatrix750[names(gse_list_main)[j], names(gse_list_main)[i]] = cor(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
  }
}


#CORPVAL
for (i in 1:length(gse_list_main)){
  for (j in i:length(gse_list_main)){
    topA = gse_list_main[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = gse_list_main[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    corpval750[names(gse_list_main)[i], names(gse_list_main)[j]] = as.numeric(cor.test(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    corpval750[names(gse_list_main)[j], names(gse_list_main)[i]] = as.numeric(cor.test(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
  }
}

#CORADJPVAL
vec = as.vector(corpval750[upper.tri(corpval750, diag = F)])
vec = p.adjust(vec, method = "BH")
tempmatrix = matrix(0, length(gse_list_main), length(gse_list_main))
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

saveRDS(logFCmatrixregr, file = "logFCmatrixregr_all.Rds")
saveRDS(SEmatrixregr, file = "SEmatrixregr_all.Rds")
saveRDS(cortestsign750, file = "cortestsign750_all.Rds")
saveRDS(totalrownamematrix, file = "totalrownamematrix_all.Rds")

# plot exapmles:
kekmatrix = cormatrix750 %>% rownames_to_column("datasetid1")
kekmatrix = gather(kekmatrix, datasetid2, corvalue, -datasetid1)
kekmatrix$corvalue = as.numeric(as.character(kekmatrix$corvalue))
kekmatrix = kekmatrix %>% filter(corvalue != 1) %>% top_n(10, corvalue) %>% distinct(corvalue, .keep_all = T)

for (i in 1:length(rownames(kekmfatrix))){
  plot1 = deming(logFCmatrix_all_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid2"]] ~ logFCmatrix_all_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid1"]] - 1)
  
  ggheatmap = ggplot(logFCmatrix_all_goodboys[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],], aes_string(x = kekmatrix[i,"datasetid1"], y = kekmatrix[i,"datasetid2"])) + geom_point() +
    geom_abline(slope = plot1$coefficients[2], intercept = 0, colour = "blue", size = 1) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_abline(slope = kres[which(colnames(logFCmatrix_all_goodboys) == kekmatrix[i, "datasetid2"])]/kres[which(colnames(logFCmatrix_all_goodboys) == kekmatrix[i, "datasetid1"])], intercept = 0, colour = "red", size = 1)
  print(ggheatmap)
  pdf(paste0("demingexample_all_", i,  ".pdf"))
  print(ggheatmap)
  dev.off()
}


#second step

deming_result_all = readRDS("deming_result_all.Rds")
logFCmatrixregr_all = readRDS("logFCmatrixregr_all.Rds")
SEmatrixregr_all = readRDS("SEmatrixregr_all.Rds")
SEmatrixregr_all = readRDS("SEmatrixregr_all.Rds")

kres = deming_result_all
for (i in 1:length(colnames(logFCmatrixregr_all))){
  SEmatrixregr_all[,i] = SEmatrixregr_all[,i] / kres[i]
  logFCmatrixregr_all[,i] = logFCmatrixregr_all[,i] / kres[i]
}
temp_colnames = colnames(logFCmatrixregr_all)
logFCmatrixregr_all = data.frame(logFCmatrixregr_all)
logFCmatrixregr_all$NACount = rowSums(is.na(logFCmatrixregr_all))
ggplot(logFCmatrixregr_all, aes(x = NACount)) + geom_density()
logFCmatrix_all_goodboys = subset(logFCmatrixregr_all, logFCmatrixregr_all$NACount <=round((2*ncol(logFCmatrixregr_all))/3))
SEmatrix_all_goodboys = subset(SEmatrixregr_all, logFCmatrixregr_all$NACount <=round((2*ncol(logFCmatrixregr_all))/3))
logFCmatrix_all_goodboys$NACount = NULL
colnames(logFCmatrix_all_goodboys) = temp_colnames


sourcedata_all = as.data.frame(colnames(logFCmatrix_all_goodboys))
rownames(sourcedata_all) = colnames(logFCmatrix_all_goodboys)
colnames(sourcedata_all) = "datasets"
#species column have made manually
sourcedata_all$species = c(rep("human", 20), rep(NA, 18))
cops = c()
for (i in strsplit(rownames(sourcedata_all), "_")){
  cops = c(cops, gsub(".Rds", "", i[length(i)]))
}
sourcedata_all[,"datasets"] = cops
all_signature = signature_builder(logFCmatrix_all_goodboys, SEmatrix_all_goodboys, sourcedata_all)
saveRDS(all_signature, file = "all_signature.Rds")

