library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
source("functions.R")
library(biomaRt)
library(preprocessCore)
library(biomaRt)
library(illuminaHumanv4.db)
library(oligoClasses)
library(org.Mm.eg.db)

####################################################################################
####################################################################################
#PLATFORM:  Illumina NextSeq 500 (Homo sapiens)
#create from RAW data, different factors and their combinations
####################################################################################
####################################################################################
current_gse = "gse127927"
current_gse_geo = getGEO(current_gse)
#STEP 0: create tables for analysis 
current_gse_phenoData = pData(current_gse_geo[[1]])
current_gse_assayData = exprs(current_gse_geo[[1]])
current_gse_featureData = fData(current_gse_geo[[1]])
current_gse_phenoData = current_gse_phenoData[6:67,]

#shitcode to deal with SRR's
pops = c()
flops = c()
for (i in list.files("Users/braaist/Documents/gse127927/")){
  print(i)
  pops = c(pops, substr(i, 1, 10))
}
pops = data.frame(pops)
tf_vec = c()
for (i in rownames(pops)){
  tf_vec = c(tf_vec, as.numeric(i) %% 2 == 0)
}
flops = pops[tf_vec,]

flops = flops[1:62]
rownames(current_gse_phenoData) = flops
current_gse_geo = Download_raw_reads("/Users/braaist/Documents/gse127927/", phenodata = current_gse_phenoData)
current_gse_phenoData = pData(current_gse_geo)
current_gse_assayData = exprs(current_gse_geo)
current_gse_featureData = fData(current_gse_geo)

#and check density
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before all", current_gse, sep = " "))


#STEP 1: logtransformation (data is not transformed)
#seems already transformed, it's actually strange
current_gse_assayData = log(current_gse_assayData + 1)

#STEP 2: subset data
#current_gse_featureData = current_gse_featureData
#rownames(current_gse_assayData) = current_gse_featureData$ID
#STEP 3: remove outliers
#current_gse_assayData = current_gse_assayData[,c(1:2, 4:10)]
#current_gse_phenoData = current_gse_phenoData[c(1:2, 4:10),]
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before", current_gse, sep = " "))


#STEP 4: filter by treshold
filter_table = current_gse_assayData >= 2
filter_table = data.frame(filter_table)
filter_table = filter_table %>% mutate(sum = rowSums(.[1:ncol(filter_table)])) %>% mutate(filter = sum>15)
filter_table = data.frame(filter_table)
filter_treshold = filter_table$filter
current_gse_assayData = current_gse_assayData[filter_treshold,]
current_gse_featureData = current_gse_featureData[filter_treshold,]
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density before", current_gse, sep = " "))

#STEP 5: scaling and quantile normalization
#scaling
current_gse_assayData = scale(current_gse_assayData)
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after scaling", current_gse, sep = " "))
colnames_ad = current_gse_phenoData$geo_accession
rownames_ad = rownames(current_gse_assayData)
#quantile normalization
current_gse_assayData = normalize.quantiles(as.matrix(current_gse_assayData))
#and check density again
current_gse_assayData_fd = data.frame(current_gse_assayData)
visualstack = stack(current_gse_assayData_fd)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("density after quantile normalization", current_gse, sep = " "))

colnames(current_gse_assayData) = colnames_ad
rownames(current_gse_assayData) = rownames_ad
current_gse_assayData = data.frame(current_gse_assayData)
#STEP 6: add entrez ids
#add entrez ids
entrez_vec = ensembl2entrez_mm(rownames(current_gse_assayData))
current_gse_assayData$entrez = entrez_vec


#delete genes without ENTREZID
NAfilter = !is.na(current_gse_assayData$entrez)
current_gse_assayData = current_gse_assayData[NAfilter,]
#STEP 7: replace genes with the same entrez id 
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_assayData)){
  print(i)
  if ((length(which(current_gse_assayData[i,'entrez'] == current_gse_assayData[,'entrez'])) != 1) && (list(which(current_gse_assayData[i,'entrez'] == current_gse_assayData[,'entrez'])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_assayData[i,'entrez'] == current_gse_assayData[,'entrez'])))
  }
}
#get sum of the RNA seq data for same ENTREZID
for (i in same_entrez_list){
  same_entrez_matrix = double()
  for (j in i){
    same_entrez_matrix = rbind(same_entrez_matrix, as.double(current_gse_assayData[j,]))}
  sum_of_samples = apply(same_entrez_matrix, 2, sum)
  current_gse_assayData[i[1],] = sum_of_samples
  i = i[-1]
  for (m in i){
    current_gse_assayData[m,] = NA
    current_gse_featureData[m,] = NA}
}

current_gse_assayData = na.omit(current_gse_assayData)
rownames(current_gse_assayData) = current_gse_assayData$entrez
current_gse_assayData$entrez = NULL
entrez_vec = rownames(current_gse_assayData)
#STEP 8: PCA plot
gse_for_cluster = scale(t(current_gse_assayData))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2])) +
  geom_point(aes(color = current_gse_phenoData$title)) + 
  labs(colour = "days", title = paste("PCA plot","", current_gse, sep = " "))
cluster_plot

#STEP 8: final scaling
#scaling
current_gse_assayData = scale(current_gse_assayData)
#and check density again
current_gse_assayData_fd = current_gse_assayData
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + 
  labs(colour = "days", title = paste("final density", current_gse, sep = " "))

current_gse_phenoData_back = current_gse_phenoData
current_gse_assayData_back = current_gse_assayData
#STEP 9: limma
#Glis1
current_gse_phenoData = current_gse_phenoData_back
current_gse_assayData = current_gse_assayData_back
current_gse_phenoData1 = current_gse_phenoData[16:20,]
current_gse_assayData1 = current_gse_assayData[,16:20]
flops = strsplit(as.character(current_gse_phenoData1$title), "-")
days_vec = c()
for (i in flops){
  days_vec = c(days_vec, i[3])
}
days_vec = sub("D", "", days_vec)
days_vec = as.numeric(days_vec)
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData1)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData1, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData1[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
colnames(first_entity_df) = c("first_entity_df", 'day')

ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("Sall4", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)
save(top_genes_limma,file="Users/braaist/Documents/top_genes_Mkk6.Rda")



#################################################################################################
#################################################################################################
#All factors IPSC
#################################################################################################
#################################################################################################
current_gse_phenoData = current_gse_phenoData_back
current_gse_assayData = current_gse_assayData_back
current_gse_phenoData1 = current_gse_phenoData[c(45:47, 53:58),]
current_gse_assayData1 = current_gse_assayData[,c(45:47, 53:58)]
days_vec = c(rep(0,3), rep(1,6))
design_matrix = cbind(c(rep(1, length(days_vec))), days_vec)
colnames(design_matrix) <- c("Intercept", "Age")
rownames(design_matrix) = colnames(current_gse_assayData1)
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
top_genes_limma = limma_maker(current_gse_assayData1, design_matrix, contrast_matrix)

first_entity = rownames(top_genes_limma[1,])
first_entity_df = current_gse_assayData1[first_entity,]
first_entity_df = data.frame(first_entity_df)
first_entity_df$day = days_vec
colnames(first_entity_df) = c("first_entity_df", 'day')

ggplot(first_entity_df, aes(day, first_entity_df, color = day))  + geom_point() +
  labs(colour = "days", x = "age", y = "expression", title = paste("Y'F", first_entity, current_gse, sep = " ")) +
  scale_x_continuous(labels = as.character(days_vec), breaks = days_vec)
save(top_genes_limma,file="Users/braaist/Documents/top_genes_YF.Rda")






