#get GSE file create phenodata table
library(GEOquery)
library(dplyr)
library(mogene10sttranscriptcluster.db)
library(ggplot2)
current_gse = getGEO("gse55272")
#create phenodata table
current_gse_phenoData = data.frame(pData(current_gse[[1]])[,c(44,41,40,42,43,25)])
current_gse_phenoData = subset(current_gse_phenoData, tissue.ch1=="Liver")
current_gse_phenoData = current_gse_phenoData %>% mutate("diet" = "-", "drug" = "-", "dose" = "-","age_intervention" = "-")
#create assayData table
current_gse_assayData = exprs(current_gse[[1]])
#get Liver tissue assayData
current_gse_assayData = current_gse_assayData[,c(1:6,19:24)]
#make norm plots(without log transformation, because current data is already logtransformed)
lognorm_current_gse_assayData = scale(current_gse_assayData)
den = apply(lognorm_current_gse_assayData, 2, density)
plot(NA, xlim=range(sapply(den, "[", "x")), ylim=range(sapply(den, "[", "y")))
mapply(lines, den, col=1:length(den))
#PCA
gse_for_cluster = scale(t(current_gse_assayData))
#check scaling
sd(gse_for_cluster[1,])
#make PCA plot
pcamodel = prcomp(gse_for_cluster)
cluster_values = as.data.frame(pcamodel$x)
cluster_plot = ggplot(cluster_values, aes(cluster_values[,1], cluster_values[,2]))
color_phd = paste(current_gse_phenoData[,2],current_gse_phenoData[,3],sep="/")
cluster_plot = cluster_plot + geom_point(aes(color = color_phd))
cluster_plot
#add ENTREZID column
x = mogene10sttranscriptclusterENTREZID
mapped_probes = mappedkeys(x)
xx = as.list(x[mapped_probes])
entrez_converter <- function(ACCNUM){
#function for convertation ACCNUM to ENTREZID
  if(is.null(xx[[as.character(ACCNUM)]])){
    return(NA)
  }else{
    return(xx[[as.character(ACCNUM)]])}
}
current_gse_featureData = fData(current_gse[[1]])
current_gse_featureData["ENTREZID"] = sapply(current_gse_featureData[,1], entrez_converter)
#delete genes without ENTREZID
NAfilter = !is.na(current_gse_featureData[,13])
current_gse_assayData = current_gse_assayData[NAfilter,]
current_gse_featureData = current_gse_featureData[NAfilter,]
#get lists of genes with the same ENTREZID
same_entrez_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(current_gse_featureData)){
  if ((length(which(current_gse_featureData[i,13] == current_gse_featureData[,13])) != 1) && (list(which(current_gse_featureData[i,13] == current_gse_featureData[,13])) %notin% same_entrez_list)){
    same_entrez_list = append(same_entrez_list, list(which(current_gse_featureData[i,13] == current_gse_featureData[,13])))
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
    current_gse_assayData[m,] = NA}
}
current_gse_assayData = na.omit(current_gse_assayData)
#draft
