library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
source("functions.R")
library(biomaRt)
library(org.Mm.eg.db)
library(preprocessCore)
library(biomaRt)

current_gse = "gse89455"
current_gse_geo = getGEO(current_gse)
#STEP 0: create phenodata table
current_gse_phenoData = pData(current_gse_geo[[1]])
current_gse_assayData = exprs(current_gse_geo[[1]])
current_gse_featureData = fData(current_gse_geo[[1]])

#STEP 0.5
#get list of genes according to agilent platform ids
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "asia")
genes = getBM(attributes=c('agilent_sureprint_g3_ge_8x60k','entrezgene_id'), mart = ensembl)


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF HDF
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: create density plot
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) +
  labs(colour = "days", title = paste("density before preprocessing", current_gse, sep = " "))

#STEP 2.5: subset data
current_gse_phenoData_HDF = current_gse_phenoData[1:21,]
current_gse_assayData_HDF = current_gse_assayData[,1:21]
current_gse_featureData_HDF = current_gse_featureData
rownames(current_gse_assayData_HDF) = current_gse_featureData_HDF$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HDF)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after removing outliers", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_HDF >= -4
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>6)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_HDF = current_gse_assayData_HDF[filter_treshold,]
current_gse_featureData_HDF = current_gse_featureData_HDF[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HDF)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold HDF", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_HDF = scale(current_gse_assayData_HDF)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HDF)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_HDF = normalize.quantiles(as.matrix(current_gse_assayData_HDF))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HDF)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization", current_gse, sep = " "))
colnames(current_gse_assayData_HDF) = rownames(current_gse_phenoData_HDF)
rownames(current_gse_assayData_HDF) = current_gse_featureData_HDF$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_HDF = data.frame(current_gse_assayData_HDF)
current_gse_assayData_HDF$entrez = NA
for (i in rownames(current_gse_assayData_HDF)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_HDF[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_HDF[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_HDF$entrez)
current_gse_assayData_HDF = current_gse_assayData_HDF[NAfilter,]
current_gse_featureData_HDF = current_gse_featureData_HDF[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_HDF)){
  print(i)
  if ((length(which(current_gse_assayData_HDF[i,'entrez'] == current_gse_assayData_HDF[,'entrez'])) != 1) && (list(which(current_gse_assayData_HDF[i,'entrez'] == current_gse_assayData_HDF[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_HDF[i,'entrez'] == current_gse_assayData_HDF[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_HDF[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_HDF[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_HDF[m,] = NA
    current_gse_featureData_HDF[m,] = NA}
}

current_gse_assayData_HDF = na.omit(current_gse_assayData_HDF)
current_gse_assayData_HDF = data.frame(current_gse_assayData_HDF)
current_gse_assayData_HDF = data.frame(current_gse_assayData_HDF)
rownames(current_gse_assayData_HDF) = current_gse_assayData_HDF$entrez
current_gse_assayData_HDF$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_HDF))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_HDF$title)) + 
  labs(colour = "days", title = paste("PCA plot","HDF", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling"
current_gse_assayData_HDF = scale(current_gse_assayData_HDF)
#and check density again
current_gse_assayData_fd = current_gse_assayData
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density", current_gse, sep = " "))


#STEP 9: limma
days_vec = strsplit(as.character(current_gse_phenoData_HDF$title), '_')
temp = c()
for (i in days_vec){
  temp = c(temp, i[2])
}
temp = sub('d', '', temp)
days_vec = as.numeric(temp)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_HDF)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_HDF, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_HDF[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID HDF", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_HDF.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF HAdMSC
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_HAdMSC = current_gse_phenoData[22:42,]
current_gse_assayData_HAdMSC = current_gse_assayData[,22:42]
current_gse_featureData_HAdMSC = current_gse_featureData
rownames(current_gse_assayData_HAdMSC) = current_gse_featureData_HAdMSC$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HAdMSC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before HAdMSC", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_HAdMSC >= -4
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>6)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_HAdMSC = current_gse_assayData_HAdMSC[filter_treshold,]
current_gse_featureData_HAdMSC = current_gse_featureData_HAdMSC[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HAdMSC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold HAdMSC", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_HAdMSC = scale(current_gse_assayData_HAdMSC)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HAdMSC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling HAdMSC", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_HAdMSC = normalize.quantiles(as.matrix(current_gse_assayData_HAdMSC))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HAdMSC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization", current_gse, sep = " "))
colnames(current_gse_assayData_HAdMSC) = rownames(current_gse_phenoData_HAdMSC)
rownames(current_gse_assayData_HAdMSC) = current_gse_featureData_HAdMSC$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_HAdMSC = data.frame(current_gse_assayData_HAdMSC)
current_gse_assayData_HAdMSC$entrez = NA
for (i in rownames(current_gse_assayData_HAdMSC)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_HAdMSC[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_HAdMSC[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_HAdMSC$entrez)
current_gse_assayData_HAdMSC = current_gse_assayData_HAdMSC[NAfilter,]
current_gse_featureData_HAdMSC = current_gse_featureData_HAdMSC[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_HAdMSC)){
  print(i)
  if ((length(which(current_gse_assayData_HAdMSC[i,'entrez'] == current_gse_assayData_HAdMSC[,'entrez'])) != 1) && (list(which(current_gse_assayData_HAdMSC[i,'entrez'] == current_gse_assayData_HAdMSC[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_HAdMSC[i,'entrez'] == current_gse_assayData_HAdMSC[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_HAdMSC[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_HAdMSC[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_HAdMSC[m,] = NA
    current_gse_featureData_HAdMSC[m,] = NA}
}

current_gse_assayData_HAdMSC = na.omit(current_gse_assayData_HAdMSC)
current_gse_assayData_HAdMSC = data.frame(current_gse_assayData_HAdMSC)
current_gse_assayData_HAdMSC = data.frame(current_gse_assayData_HAdMSC)
rownames(current_gse_assayData_HAdMSC) = current_gse_assayData_HAdMSC$entrez
current_gse_assayData_HAdMSC$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_HAdMSC))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_HAdMSC$title)) + 
  labs(colour = "days", title = paste("PCA plot","HAdMSC", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_HAdMSC = scale(current_gse_assayData_HAdMSC)
#and check density again
current_gse_assayData_fd = current_gse_assayData_HAdMSC
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density HAdMSC", current_gse, sep = " "))


#STEP 9: limma
days_vec = strsplit(as.character(current_gse_phenoData_HAdMSC$title), '_')
temp = c()
for (i in days_vec){
  temp = c(temp, i[2])
}
temp = sub('d', '', temp)
days_vec = as.numeric(temp)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_HAdMSC)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_HAdMSC, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_HAdMSC[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID HAdMSC", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_HAdMSC.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF HA
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_HA = current_gse_phenoData[43:63,]
current_gse_assayData_HA = current_gse_assayData[,43:63]
current_gse_featureData_HA = current_gse_featureData
rownames(current_gse_assayData_HA) = current_gse_featureData_HA$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HA)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before HA", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_HA >= -4
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>6)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_HA = current_gse_assayData_HA[filter_treshold,]
current_gse_featureData_HA = current_gse_featureData_HA[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HA)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold HA", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_HA = scale(current_gse_assayData_HA)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HA)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling HA", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_HA = normalize.quantiles(as.matrix(current_gse_assayData_HA))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HA)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization", current_gse, sep = " "))
colnames(current_gse_assayData_HA) = rownames(current_gse_phenoData_HA)
rownames(current_gse_assayData_HA) = current_gse_featureData_HA$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_HA = data.frame(current_gse_assayData_HA)
current_gse_assayData_HA$entrez = NA
for (i in rownames(current_gse_assayData_HA)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_HA[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_HA[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_HA$entrez)
current_gse_assayData_HA = current_gse_assayData_HA[NAfilter,]
current_gse_featureData_HA = current_gse_featureData_HA[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_HA)){
  print(i)
  if ((length(which(current_gse_assayData_HA[i,'entrez'] == current_gse_assayData_HA[,'entrez'])) != 1) && (list(which(current_gse_assayData_HA[i,'entrez'] == current_gse_assayData_HA[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_HA[i,'entrez'] == current_gse_assayData_HA[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_HA[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_HA[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_HA[m,] = NA
    current_gse_featureData_HA[m,] = NA}
}

current_gse_assayData_HA = na.omit(current_gse_assayData_HA)
current_gse_assayData_HA = data.frame(current_gse_assayData_HA)
current_gse_assayData_HA = data.frame(current_gse_assayData_HA)
rownames(current_gse_assayData_HA) = current_gse_assayData_HA$entrez
current_gse_assayData_HA$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_HA))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_HA$title)) + 
  labs(colour = "days", title = paste("PCA plot","HA", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_HA = scale(current_gse_assayData_HA)
#and check density again
current_gse_assayData_fd = current_gse_assayData_HA
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density HA", current_gse, sep = " "))


#STEP 9: limma
days_vec = strsplit(as.character(current_gse_phenoData_HA$title), '_')
temp = c()
for (i in days_vec){
  temp = c(temp, i[2])
}
temp = sub('d', '', temp)
days_vec = as.numeric(temp)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_HA)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_HA, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_HA[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID HA", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_HA.Rda")




#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF HBEC
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_HBEC = current_gse_phenoData[64:84,]
current_gse_assayData_HBEC = current_gse_assayData[,64:84]
current_gse_featureData_HBEC = current_gse_featureData
rownames(current_gse_assayData_HBEC) = current_gse_featureData_HBEC$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HBEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before HBEC", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_HBEC >= -4
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>6)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_HBEC = current_gse_assayData_HBEC[filter_treshold,]
current_gse_featureData_HBEC = current_gse_featureData_HBEC[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HBEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold HBEC", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_HBEC = scale(current_gse_assayData_HBEC)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HBEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling HBEC", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_HBEC = normalize.quantiles(as.matrix(current_gse_assayData_HBEC))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HBEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization HBEC", current_gse, sep = " "))
colnames(current_gse_assayData_HBEC) = rownames(current_gse_phenoData_HBEC)
rownames(current_gse_assayData_HBEC) = current_gse_featureData_HBEC$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_HBEC = data.frame(current_gse_assayData_HBEC)
current_gse_assayData_HBEC$entrez = NA
for (i in rownames(current_gse_assayData_HBEC)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_HBEC[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_HBEC[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_HBEC$entrez)
current_gse_assayData_HBEC = current_gse_assayData_HBEC[NAfilter,]
current_gse_featureData_HBEC = current_gse_featureData_HBEC[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_HBEC)){
  print(i)
  if ((length(which(current_gse_assayData_HBEC[i,'entrez'] == current_gse_assayData_HBEC[,'entrez'])) != 1) && (list(which(current_gse_assayData_HBEC[i,'entrez'] == current_gse_assayData_HBEC[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_HBEC[i,'entrez'] == current_gse_assayData_HBEC[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_HBEC[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_HBEC[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_HBEC[m,] = NA
    current_gse_featureData_HBEC[m,] = NA}
}

current_gse_assayData_HBEC = na.omit(current_gse_assayData_HBEC)
current_gse_assayData_HBEC = data.frame(current_gse_assayData_HBEC)
current_gse_assayData_HBEC = data.frame(current_gse_assayData_HBEC)
rownames(current_gse_assayData_HBEC) = current_gse_assayData_HBEC$entrez
current_gse_assayData_HBEC$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_HBEC))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_HBEC$title)) + 
  labs(colour = "days", title = paste("PCA plot","HBEC", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_HBEC = scale(current_gse_assayData_HBEC)
#and check density again
current_gse_assayData_fd = current_gse_assayData_HBEC
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density HBEC", current_gse, sep = " "))


#STEP 9: limma
days_vec = strsplit(as.character(current_gse_phenoData_HBEC$title), '_')
temp = c()
for (i in days_vec){
  temp = c(temp, i[2])
}
temp = sub('d', '', temp)
days_vec = as.numeric(temp)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_HBEC)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_HBEC, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_HBEC[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID HBEC", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_HBEC.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF HPrEC
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_HPrEC = current_gse_phenoData[85:105,]
current_gse_assayData_HPrEC = current_gse_assayData[,85:105]
current_gse_featureData_HPrEC = current_gse_featureData
rownames(current_gse_assayData_HPrEC) = current_gse_featureData_HPrEC$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HPrEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before HPrEC", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_HPrEC >= -4
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>6)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_HPrEC = current_gse_assayData_HPrEC[filter_treshold,]
current_gse_featureData_HPrEC = current_gse_featureData_HPrEC[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HPrEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold HPrEC", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_HPrEC = scale(current_gse_assayData_HPrEC)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HPrEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling HPrEC", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_HPrEC = normalize.quantiles(as.matrix(current_gse_assayData_HPrEC))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_HPrEC)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization HPrEC", current_gse, sep = " "))
colnames(current_gse_assayData_HPrEC) = rownames(current_gse_phenoData_HPrEC)
rownames(current_gse_assayData_HPrEC) = current_gse_featureData_HPrEC$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_HPrEC = data.frame(current_gse_assayData_HPrEC)
current_gse_assayData_HPrEC$entrez = NA
for (i in rownames(current_gse_assayData_HPrEC)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_HPrEC[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_HPrEC[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_HPrEC$entrez)
current_gse_assayData_HPrEC = current_gse_assayData_HPrEC[NAfilter,]
current_gse_featureData_HPrEC = current_gse_featureData_HPrEC[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_HPrEC)){
  print(i)
  if ((length(which(current_gse_assayData_HPrEC[i,'entrez'] == current_gse_assayData_HPrEC[,'entrez'])) != 1) && (list(which(current_gse_assayData_HPrEC[i,'entrez'] == current_gse_assayData_HPrEC[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_HPrEC[i,'entrez'] == current_gse_assayData_HPrEC[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_HPrEC[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_HPrEC[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_HPrEC[m,] = NA
    current_gse_featureData_HPrEC[m,] = NA}
}

current_gse_assayData_HPrEC = na.omit(current_gse_assayData_HPrEC)
current_gse_assayData_HPrEC = data.frame(current_gse_assayData_HPrEC)
current_gse_assayData_HPrEC = data.frame(current_gse_assayData_HPrEC)
rownames(current_gse_assayData_HPrEC) = current_gse_assayData_HPrEC$entrez
current_gse_assayData_HPrEC$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_HPrEC))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_HPrEC$title)) + 
  labs(colour = "days", title = paste("PCA plot","HPrEC", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_HPrEC = scale(current_gse_assayData_HPrEC)
#and check density again
current_gse_assayData_fd = current_gse_assayData_HPrEC
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density HPrEC", current_gse, sep = " "))


#STEP 9: limma
days_vec = strsplit(as.character(current_gse_phenoData_HPrEC$title), '_')
temp = c()
for (i in days_vec){
  temp = c(temp, i[2])
}
temp = sub('d', '', temp)
days_vec = as.numeric(temp)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_HPrEC)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_HPrEC, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_HPrEC[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID HPrEC", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_HPrEC.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF O versus Mock
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_O = current_gse_phenoData[c(106:108,151:153),]
current_gse_assayData_O = current_gse_assayData[,c(106:108,151:153)]
current_gse_featureData_O = current_gse_featureData
rownames(current_gse_assayData_O) = current_gse_featureData_O$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_O)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before O", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_O >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_O = current_gse_assayData_O[filter_treshold,]
current_gse_featureData_O = current_gse_featureData_O[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_O)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold O", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_O = scale(current_gse_assayData_O)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_O)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling O", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_O = normalize.quantiles(as.matrix(current_gse_assayData_O))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_O)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization O", current_gse, sep = " "))
colnames(current_gse_assayData_O) = rownames(current_gse_phenoData_O)
rownames(current_gse_assayData_O) = current_gse_featureData_O$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_O = data.frame(current_gse_assayData_O)
current_gse_assayData_O$entrez = NA
for (i in rownames(current_gse_assayData_O)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_O[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_O[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_O$entrez)
current_gse_assayData_O = current_gse_assayData_O[NAfilter,]
current_gse_featureData_O = current_gse_featureData_O[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_O)){
  print(i)
  if ((length(which(current_gse_assayData_O[i,'entrez'] == current_gse_assayData_O[,'entrez'])) != 1) && (list(which(current_gse_assayData_O[i,'entrez'] == current_gse_assayData_O[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_O[i,'entrez'] == current_gse_assayData_O[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_O[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_O[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_O[m,] = NA
    current_gse_featureData_O[m,] = NA}
}

current_gse_assayData_O = na.omit(current_gse_assayData_O)
current_gse_assayData_O = data.frame(current_gse_assayData_O)
current_gse_assayData_O = data.frame(current_gse_assayData_O)
rownames(current_gse_assayData_O) = current_gse_assayData_O$entrez
current_gse_assayData_O$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_O))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_O$title)) + 
  labs(colour = "days", title = paste("PCA plot","O", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_O = scale(current_gse_assayData_O)
#and check density again
current_gse_assayData_fd = current_gse_assayData_O
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density O", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_O)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_O, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_O[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID O", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_O.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF S versus Mock
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_S = current_gse_phenoData[c(109:111,151:153),]
current_gse_assayData_S = current_gse_assayData[,c(109:111,151:153)]
current_gse_featureData_S = current_gse_featureData
rownames(current_gse_assayData_S) = current_gse_featureData_S$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_S)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before S", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_S >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_S = current_gse_assayData_S[filter_treshold,]
current_gse_featureData_S = current_gse_featureData_S[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_S)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold S", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_S = scale(current_gse_assayData_S)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_S)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling O", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_S = normalize.quantiles(as.matrix(current_gse_assayData_S))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_S)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization O", current_gse, sep = " "))
colnames(current_gse_assayData_S) = rownames(current_gse_phenoData_S)
rownames(current_gse_assayData_S) = current_gse_featureData_S$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_S = data.frame(current_gse_assayData_S)
current_gse_assayData_S$entrez = NA
for (i in rownames(current_gse_assayData_S)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_S[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_S[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_S$entrez)
current_gse_assayData_S = current_gse_assayData_S[NAfilter,]
current_gse_featureData_S = current_gse_featureData_S[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_S)){
  print(i)
  if ((length(which(current_gse_assayData_S[i,'entrez'] == current_gse_assayData_S[,'entrez'])) != 1) && (list(which(current_gse_assayData_S[i,'entrez'] == current_gse_assayData_S[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_S[i,'entrez'] == current_gse_assayData_S[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_S[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_S[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_S[m,] = NA
    current_gse_featureData_S[m,] = NA}
}

current_gse_assayData_S = na.omit(current_gse_assayData_S)
current_gse_assayData_S = data.frame(current_gse_assayData_S)
current_gse_assayData_S = data.frame(current_gse_assayData_S)
rownames(current_gse_assayData_S) = current_gse_assayData_S$entrez
current_gse_assayData_S$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_S))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_S$title)) + 
  labs(colour = "days", title = paste("PCA plot","S", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_S = scale(current_gse_assayData_S)
#and check density again
current_gse_assayData_fd = current_gse_assayData_S
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density S", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_S)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_S, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_S[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID S", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_S.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF K versus Mock
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_K = current_gse_phenoData[c(112:114,151:153),]
current_gse_assayData_K = current_gse_assayData[,c(112:114,151:153)]
current_gse_featureData_K = current_gse_featureData
rownames(current_gse_assayData_K) = current_gse_featureData_K$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_K)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before K", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_K >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_K = current_gse_assayData_K[filter_treshold,]
current_gse_featureData_K = current_gse_featureData_K[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_K)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold K", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_K = scale(current_gse_assayData_K)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_K)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling K", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_K = normalize.quantiles(as.matrix(current_gse_assayData_K))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_K)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization K", current_gse, sep = " "))
colnames(current_gse_assayData_K) = rownames(current_gse_phenoData_K)
rownames(current_gse_assayData_K) = current_gse_featureData_K$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_K = data.frame(current_gse_assayData_K)
current_gse_assayData_K$entrez = NA
for (i in rownames(current_gse_assayData_K)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_K[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_K[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_K$entrez)
current_gse_assayData_K = current_gse_assayData_K[NAfilter,]
current_gse_featureData_K = current_gse_featureData_K[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_K)){
  print(i)
  if ((length(which(current_gse_assayData_K[i,'entrez'] == current_gse_assayData_K[,'entrez'])) != 1) && (list(which(current_gse_assayData_K[i,'entrez'] == current_gse_assayData_K[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_K[i,'entrez'] == current_gse_assayData_K[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_K[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_K[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_K[m,] = NA
    current_gse_featureData_K[m,] = NA}
}

current_gse_assayData_K = na.omit(current_gse_assayData_K)
current_gse_assayData_K = data.frame(current_gse_assayData_K)
current_gse_assayData_K = data.frame(current_gse_assayData_K)
rownames(current_gse_assayData_K) = current_gse_assayData_K$entrez
current_gse_assayData_K$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_K))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_K$title)) + 
  labs(colour = "days", title = paste("PCA plot","K", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_K = scale(current_gse_assayData_K)
#and check density again
current_gse_assayData_fd = current_gse_assayData_K
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density K", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_K)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_K, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_K[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID K", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_K.Rda")

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF M versus Mock
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_M = current_gse_phenoData[c(115:117,151:153),]
current_gse_assayData_M = current_gse_assayData[,c(115:117,151:153)]
current_gse_featureData_M = current_gse_featureData
rownames(current_gse_assayData_M) = current_gse_featureData_M$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_M)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before M", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_M >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_M = current_gse_assayData_M[filter_treshold,]
current_gse_featureData_M = current_gse_featureData_M[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_M)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold M", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_M = scale(current_gse_assayData_M)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_M)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling M", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_M = normalize.quantiles(as.matrix(current_gse_assayData_M))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_M)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization M", current_gse, sep = " "))
colnames(current_gse_assayData_M) = rownames(current_gse_phenoData_M)
rownames(current_gse_assayData_M) = current_gse_featureData_M$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_M = data.frame(current_gse_assayData_M)
current_gse_assayData_M$entrez = NA
for (i in rownames(current_gse_assayData_M)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_M[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_M[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_M$entrez)
current_gse_assayData_M = current_gse_assayData_M[NAfilter,]
current_gse_featureData_M = current_gse_featureData_M[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_M)){
  print(i)
  if ((length(which(current_gse_assayData_M[i,'entrez'] == current_gse_assayData_M[,'entrez'])) != 1) && (list(which(current_gse_assayData_M[i,'entrez'] == current_gse_assayData_M[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_M[i,'entrez'] == current_gse_assayData_M[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_M[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_M[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_M[m,] = NA
    current_gse_featureData_M[m,] = NA}
}

current_gse_assayData_M = na.omit(current_gse_assayData_M)
current_gse_assayData_M = data.frame(current_gse_assayData_M)
current_gse_assayData_M = data.frame(current_gse_assayData_M)
rownames(current_gse_assayData_M) = current_gse_assayData_M$entrez
current_gse_assayData_M$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_M))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_M$title)) + 
  labs(colour = "days", title = paste("PCA plot","M", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_M = scale(current_gse_assayData_M)
#and check density again
current_gse_assayData_fd = current_gse_assayData_M
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density M", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_M)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_M, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_M[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID M", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_M.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF OS versus Mock
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_OS = current_gse_phenoData[c(118:120,151:153),]
current_gse_assayData_OS = current_gse_assayData[,c(118:120,151:153)]
current_gse_featureData_OS = current_gse_featureData
rownames(current_gse_assayData_OS) = current_gse_featureData_OS$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OS)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before OS", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_OS >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_OS = current_gse_assayData_OS[filter_treshold,]
current_gse_featureData_OS = current_gse_featureData_OS[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OS)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold OS", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_OS = scale(current_gse_assayData_OS)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OS)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling OS", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_OS = normalize.quantiles(as.matrix(current_gse_assayData_OS))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OS)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization OS", current_gse, sep = " "))
colnames(current_gse_assayData_OS) = rownames(current_gse_phenoData_OS)
rownames(current_gse_assayData_OS) = current_gse_featureData_OS$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_OS = data.frame(current_gse_assayData_OS)
current_gse_assayData_OS$entrez = NA
for (i in rownames(current_gse_assayData_OS)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_OS[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_OS[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_OS$entrez)
current_gse_assayData_OS = current_gse_assayData_OS[NAfilter,]
current_gse_featureData_OS = current_gse_featureData_OS[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_OS)){
  print(i)
  if ((length(which(current_gse_assayData_OS[i,'entrez'] == current_gse_assayData_OS[,'entrez'])) != 1) && (list(which(current_gse_assayData_OS[i,'entrez'] == current_gse_assayData_OS[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_OS[i,'entrez'] == current_gse_assayData_OS[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_OS[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_OS[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_OS[m,] = NA
    current_gse_featureData_OS[m,] = NA}
}

current_gse_assayData_OS = na.omit(current_gse_assayData_OS)
current_gse_assayData_OS = data.frame(current_gse_assayData_OS)
current_gse_assayData_OS = data.frame(current_gse_assayData_OS)
rownames(current_gse_assayData_OS) = current_gse_assayData_OS$entrez
current_gse_assayData_OS$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_OS))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_OS$title)) + 
  labs(colour = "days", title = paste("PCA plot","OS", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_OS = scale(current_gse_assayData_OS)
#and check density again
current_gse_assayData_fd = current_gse_assayData_OS
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density OS", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_OS)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_OS, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_OS[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID OS", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_OS.Rda")



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF OK versus Mock
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_OK = current_gse_phenoData[c(121:123,151:153),]
current_gse_assayData_OK = current_gse_assayData[,c(121:123,151:153)]
current_gse_featureData_OK = current_gse_featureData
rownames(current_gse_assayData_OK) = current_gse_featureData_OK$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OK)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before OK", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_OK >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_OK = current_gse_assayData_OK[filter_treshold,]
current_gse_featureData_OK = current_gse_featureData_OK[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OK)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold OK", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_OK = scale(current_gse_assayData_OK)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OK)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling OK", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_OK = normalize.quantiles(as.matrix(current_gse_assayData_OK))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OK)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization OK", current_gse, sep = " "))
colnames(current_gse_assayData_OK) = rownames(current_gse_phenoData_OK)
rownames(current_gse_assayData_OK) = current_gse_featureData_OK$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_OK = data.frame(current_gse_assayData_OK)
current_gse_assayData_OK$entrez = NA
for (i in rownames(current_gse_assayData_OK)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_OK[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_OK[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_OK$entrez)
current_gse_assayData_OK = current_gse_assayData_OK[NAfilter,]
current_gse_featureData_OK = current_gse_featureData_OK[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_OK)){
  print(i)
  if ((length(which(current_gse_assayData_OK[i,'entrez'] == current_gse_assayData_OK[,'entrez'])) != 1) && (list(which(current_gse_assayData_OK[i,'entrez'] == current_gse_assayData_OK[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_OK[i,'entrez'] == current_gse_assayData_OK[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_OK[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_OK[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_OK[m,] = NA
    current_gse_featureData_OK[m,] = NA}
}

current_gse_assayData_OK = na.omit(current_gse_assayData_OK)
current_gse_assayData_OK = data.frame(current_gse_assayData_OK)
current_gse_assayData_OK = data.frame(current_gse_assayData_OK)
rownames(current_gse_assayData_OK) = current_gse_assayData_OK$entrez
current_gse_assayData_OK$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_OK))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_OK$title)) + 
  labs(colour = "days", title = paste("PCA plot","OK", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_OK = scale(current_gse_assayData_OK)
#and check density again
current_gse_assayData_fd = current_gse_assayData_OK
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density OK", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_OK)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_OK, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_OK[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID OK", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_OK.Rda")

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF others factors is the same
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


#STEP 1: logtransformation (data is not transformed)
#current_gse_assayData = log(current_gse_assayData+ 1)

#STEP 2: subset data
current_gse_phenoData_OSKM = current_gse_phenoData[c(148:150,151:153),]
current_gse_assayData_OSKM = current_gse_assayData[,c(148:150,151:153)]
current_gse_featureData_OSKM = current_gse_featureData
rownames(current_gse_assayData_OSKM) = current_gse_featureData_OSKM$ProbeName
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OSKM)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before OSKM", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData_OSKM >= -5
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData_OSKM = current_gse_assayData_OSKM[filter_treshold,]
current_gse_featureData_OSKM = current_gse_featureData_OSKM[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OSKM)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold OSKM", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData_OSKM = scale(current_gse_assayData_OSKM)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OSKM)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling OSKM", current_gse, sep = " "))

#quantile normalization
current_gse_assayData_OSKM = normalize.quantiles(as.matrix(current_gse_assayData_OSKM))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData_OSKM)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization OSKM", current_gse, sep = " "))
colnames(current_gse_assayData_OSKM) = rownames(current_gse_phenoData_OSKM)
rownames(current_gse_assayData_OSKM) = current_gse_featureData_OSKM$ProbeName

#STEP 6: add entrez ids
#add entrez ids
current_gse_assayData_OSKM = data.frame(current_gse_assayData_OSKM)
current_gse_assayData_OSKM$entrez = NA
for (i in rownames(current_gse_assayData_OSKM)){
  print(i)
  if (i %in% genes$agilent_sureprint_g3_ge_8x60k &&
      length(genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]) == 1){
    current_gse_assayData_OSKM[i, 'entrez'] = genes[genes$agilent_sureprint_g3_ge_8x60k == i, 2]
  }else{
    current_gse_assayData_OSKM[i, 'entrez'] = NA
  }
}

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData_OSKM$entrez)
current_gse_assayData_OSKM = current_gse_assayData_OSKM[NAfilter,]
current_gse_featureData_OSKM = current_gse_featureData_OSKM[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData_OSKM)){
  print(i)
  if ((length(which(current_gse_assayData_OSKM[i,'entrez'] == current_gse_assayData_OSKM[,'entrez'])) != 1) && (list(which(current_gse_assayData_OSKM[i,'entrez'] == current_gse_assayData_OSKM[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData_OSKM[i,'entrez'] == current_gse_assayData_OSKM[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData_OSKM[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData_OSKM[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData_OSKM[m,] = NA
    current_gse_featureData_OSKM[m,] = NA}
}

current_gse_assayData_OSKM = na.omit(current_gse_assayData_OSKM)
current_gse_assayData_OSKM = data.frame(current_gse_assayData_OSKM)
current_gse_assayData_OSKM = data.frame(current_gse_assayData_OSKM)
rownames(current_gse_assayData_OSKM) = current_gse_assayData_OSKM$entrez
current_gse_assayData_OSKM$entrez = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData_OSKM))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData_OSKM$title)) + 
  labs(colour = "days", title = paste("PCA plot","OSKM", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData_OSKM = scale(current_gse_assayData_OSKM)
#and check density again
current_gse_assayData_fd = current_gse_assayData_OSKM
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density OSKM", current_gse, sep = " "))


#STEP 9: limma
days_vec = c(rep(1,3), rep(0,3))
temp = sub('d', '', temp)
design_matrix = cbind(c(rep(1, 6)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData_OSKM)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData_OSKM, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData_OSKM[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID OSKM", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse89455/top_genes_OSKM.Rda")








