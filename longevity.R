#get GSE file create phenodata table
library(GEOquery)
library(dplyr)
library(mogene10sttranscriptcluster.db)
current_gse = getGEO("gse55272")
#create phenodata table
current_gse_phenodata = data.frame(pData(phenoData(current_gse[[1]]))[,c(44,41,40,42,43,25)])
current_gse_phenodata = subset(current_gse_phenodata, tissue.ch1=="Liver")
current_gse_phenodata = current_gse_phenodata %>% mutate("diet" = "-", "drug" = "-", "dose" = "-","age_intervention" = "-")
#create assayData table
current_gse_assayData = assayData((current_gse[[1]]))[["exprs"]]
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
color_phd = paste(current_gse_phenodata[,2],current_gse_phenodata[,3],sep="/")
cluster_plot = cluster_plot + geom_point(aes(color = color_phd))
cluster_plot
#add ENTREZID column

x = mogene10sttranscriptclusterENTREZID
mapped_probes = mappedkeys(x)
xx = as.list(x[mapped_probes])

entrez_converter <- function(ACCNUM){
#function for convertation ACCNUM to ENTREZID
  + if(is.null(xx[[as.character(ACCNUM)]])){
    + return(NA)}
  + else{
    + return(xx[[as.character(ACCNUM)]])}}
full_fdata["ENTREZID"] = sapply(full_fdata[,1], entrez_converter)


