library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
source("functions.R")
library(biomaRt)
library(preprocessCore)
library(biomaRt)
library(illuminaHumanv3.db)
####################################################################################
####################################################################################
#PLATFORM:  Illumina HumanRef-8 v3.0 expression beadchip
#HFF
####################################################################################
####################################################################################


current_gse = "gse28688"
current_gse_geo = getGEO(current_gse)
#STEP 0: create phenodata table
current_gse_phenoData = pData(current_gse_geo[[1]])
current_gse_assayData = exprs(current_gse_geo[[1]])
current_gse_featureData = fData(current_gse_geo[[1]])
#and check density
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before all", current_gse, sep = " "))



################################################################################################
################################################################################################
#HFF
################################################################################################
################################################################################################


#STEP 1: logtransformation (data is not transformed)
current_gse_assayData = log(current_gse_assayData + 1)

#STEP 2: subset data
current_gse_phenoData = current_gse_phenoData[1:8,]
current_gse_assayData = current_gse_assayData[,1:8]
current_gse_featureData = current_gse_featureData
rownames(current_gse_assayData) = current_gse_featureData$ID
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,c(1:2, 4:10)]
#current_gse_phenoData = current_gse_phenoData[c(1:2, 4:10),]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before HFF", current_gse, sep = " "))



#STEP 4: filter by treshold
#looks quite good, no filtering
#filter_table = current_gse_assayData_MEF >= -5
#filter_table = data.frame(filter_table)
#filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>2)
#filter_table = data.frame(filter_table)
#filter_treshold = filter_table$filter
#current_gse_assayData_MEF = current_gse_assayData_MEF[filter_treshold,]
#current_gse_featureData_MEF = current_gse_featureData_MEF[filter_treshold,]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after filtering by treshold ", current_gse, sep = " "))


#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData = scale(current_gse_assayData)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling HFF", current_gse, sep = " "))

#quantile normalization
current_gse_assayData = normalize.quantiles(as.matrix(current_gse_assayData))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization HFF", current_gse, sep = " "))
colnames(current_gse_assayData) = rownames(current_gse_phenoData)
rownames(current_gse_assayData) = current_gse_featureData$ID

#STEP 6: add entrez ids
#add entrez ids
entrez_vec = ilmnav32entrez_hs(current_gse_featureData$ID)
current_gse_featureData$entrez = entrez_vec

#delete genes without ENTREZID
NAfilter = !is.na(current_gse_featureData$entrez)
current_gse_assayData = current_gse_assayData[NAfilter,]
current_gse_featureData = current_gse_featureData[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_featureData)){
  print(i)
  if ((length(which(current_gse_featureData[i,'entrez'] == current_gse_featureData[,'entrez'])) != 1) && (list(which(current_gse_featureData[i,'entrez'] == current_gse_featureData[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_featureData[i,'entrez'] == current_gse_featureData[,'entrez'])))
  }
}
#get mean of the microarray data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData[j,]))}
  mean_of_samples = apply(same_entrez_matrix, 2, mean)
  current_gse_assayData[i[1],] = mean_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData[m,] = NA
    current_gse_featureData[m,] = NA}
}

current_gse_assayData = na.omit(current_gse_assayData)
current_gse_featureData = current_gse_featureData[!is.na(current_gse_featureData$entrez), ]
current_gse_assayData = data.frame(current_gse_assayData)
rownames(current_gse_assayData) = current_gse_featureData$entrez

#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData$title)) + 
  labs(colour = "days", title = paste("PCA plot","HFF cells", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData = scale(current_gse_assayData)
#and check density again
current_gse_assayData_fd = current_gse_assayData
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density HFF", current_gse, sep = " "))
colnames(current_gse_assayData) = current_gse_phenoData$geo_accession

#STEP 9: limma

days_vec = c(0,0,1,1,2,2,3,3)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("entrez ID HFF", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)

save(top_genes_limma,file="gse28688/top_genes_HFF.Rda")



