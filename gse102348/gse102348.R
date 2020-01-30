library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
source("functions.R")
library(biomaRt)
library(org.Mm.eg.db)
library(preprocessCore)

current_gse = "gse102348"
current_gse_geo = getGEO(current_gse)
#STEP 0: create phenodata table
current_gse_phenoData = pData(current_gse_geo[[1]])
current_gse_assayData = exprs(current_gse_geo[[1]])
current_gse_featureData = fData(current_gse_geo[[1]])

#There is no available assayData and featureData in this dataset, so we need to get it from site
#Only WT were used for analysis


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF FIRST WT FILE
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
current_gse_assayData = read.csv("Documents/longevity/longevity/gse102348/WT-1_fpkm.txt", sep = '\t')


#STEP 1: logtransformation (data is not transformed)
current_gse_assayData[,3:9] = log(current_gse_assayData[,3:9] + 1)

#STEP 2: create density plot
current_gse_assayData_fd = current_gse_assayData[,3:9]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) +
  labs(colour = "days", title = paste("density before preprocessing", current_gse, sep = " "))


#STEP 3: remove outliers
current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after removing outliers", current_gse, sep = " "))
  


#STEP 4: filter by treshold
filter_table = current_gse_assayData[,3:7] >= 0.2
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData = current_gse_assayData[filter_treshold,]
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData[,3:7] = scale(current_gse_assayData[,3:7])
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling", current_gse, sep = " "))

#quantile normalization
current_gse_assayData[,3:7] = normalize.quantiles(as.matrix(current_gse_assayData[,3:7]))
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization", current_gse, sep = " "))


#STEP 6: add entrez ids
#add entrez ids
ensembl_ids = current_gse_assayData$Gene.ID
entrez_ids = ensembl2entrez_mm(ensembl_ids)
current_gse_assayData$entrez_id = entrez_ids
#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData$entrez_id)
current_gse_assayData = current_gse_assayData[NAfilter,]

#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData)){
  print(i)
  if ((length(which(current_gse_assayData[i,8] == current_gse_assayData[,8])) != 1) && (list(which(current_gse_assayData[i,8] == current_gse_assayData[,8])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData[i,8] == current_gse_assayData[,8])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
  same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData[j,3:7]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData[i[1],3:7] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData[m,] = NA}
}

current_gse_assayData = na.omit(current_gse_assayData)
current_gse_assayData = data.frame(current_gse_assayData)

#STEP 7.5 due to special features of thid dataset
rownames(current_gse_assayData) = current_gse_assayData$entrez_id
current_gse_assayData$Gene.ID = NULL
current_gse_assayData$Symbol = NULL
current_gse_assayData$entrez_id = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = colnames(current_gse_assayData))) + 
  labs(colour = "days", title = paste("PCA plot", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling"
current_gse_assayData = scale(current_gse_assayData)
#and check density again
current_gse_assayData_fd = current_gse_assayData
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density", current_gse, sep = " "))



#STEP 9: limma
days_vec = c(0, 2, 4, 6, 8)
design_matrix = cbind(c(rep(1, 5)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = rownames(first_entity_df)
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID", first_entity, current_gse, sep = " "))

save(top_genes_limma,file="Documents/longevity/longevity/gse102348//top_genes_WT1.Rda")


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#ANALYSIS OF SECOND WT FILE
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
current_gse_assayData = read.csv("Documents/longevity/longevity/gse102348/WT-2_fpkm.txt", sep = '\t')

#STEP 1: logtransformation (data is not transformed)
current_gse_assayData[,3:8] = log(current_gse_assayData[,3:8] + 1)

#STEP 2: create density plot
current_gse_assayData_fd = current_gse_assayData[,3:8]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) +
  labs(colour = "days", title = paste("density before preprocessing", current_gse, sep = " "))


#STEP 3: remove outliers
current_gse_assayData = current_gse_assayData[,1:7]
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after removing outliers", current_gse, sep = " "))



#STEP 4: filter by treshold
filter_table = current_gse_assayData[,3:7] >= 0.2
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData = current_gse_assayData[filter_treshold,]
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData[,3:7] = scale(current_gse_assayData[,3:7])
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling", current_gse, sep = " "))

#quantile normalization
current_gse_assayData[,3:7] = normalize.quantiles(as.matrix(current_gse_assayData[,3:7]))
#and check density again
current_gse_assayData_fd = current_gse_assayData[,3:7]
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization", current_gse, sep = " "))


#STEP 6: add entrez ids
#add entrez ids
ensembl_ids = current_gse_assayData$Gene.ID
entrez_ids = ensembl2entrez_mm(ensembl_ids)
current_gse_assayData$entrez_id = entrez_ids
#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData$entrez_id)
current_gse_assayData = current_gse_assayData[NAfilter,]

#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData)){
  print(i)
  if ((length(which(current_gse_assayData[i,8] == current_gse_assayData[,8])) != 1) && (list(which(current_gse_assayData[i,8] == current_gse_assayData[,8])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData[i,8] == current_gse_assayData[,8])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData[j,3:7]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData[i[1],3:7] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData[m,] = NA}
}

current_gse_assayData = na.omit(current_gse_assayData)
current_gse_assayData = data.frame(current_gse_assayData)

#STEP 7.5 due to special features of thid dataset
rownames(current_gse_assayData) = current_gse_assayData$entrez_id
current_gse_assayData$Gene.ID = NULL
current_gse_assayData$Symbol = NULL
current_gse_assayData$entrez_id = NULL

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = colnames(current_gse_assayData))) + 
  labs(colour = "days", title = paste("PCA plot", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling"
current_gse_assayData = scale(current_gse_assayData)
#and check density again
current_gse_assayData_fd = current_gse_assayData
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density", current_gse, sep = " "))



#STEP 9: limma
days_vec = c(0, 2, 4, 6, 8)
design_matrix = cbind(c(rep(1, 5)), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = rownames(first_entity_df)
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID", first_entity, current_gse, sep = " "))

save(top_genes_limma,file="Documents/longevity/longevity/gse102348/top_genes_WT2.Rda")





