#get GSE file create phenodata table
library(GEOquery)
library(dplyr)
library(mogene10sttranscriptcluster.db)
library(ggplot2)
library(limma)
library(VennDiagram)
source("functions.R")
current_gse = getGEO("gse55272")
#create phenodata table
current_gse_phenoData = data.frame(pData(current_gse[[1]])[,c(44,41,40,42,43,25)])
current_gse_phenoData = subset(current_gse_phenoData, tissue.ch1=="Liver")
current_gse_phenoData = current_gse_phenoData %>% mutate("diet" = "-", "drug" = "-", "dose" = "-","age_intervention" = "-")
#create assayData table
current_gse_assayData = exprs(current_gse[[1]])
#get Liver tissue assayData
current_gse_assayData = current_gse_assayData[,c(1:6,19:24)]
#add ENTREZID column
x = mogene10sttranscriptclusterENTREZID
mapped_probes = mappedkeys(x)
xx = as.list(x[mapped_probes])
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
    current_gse_assayData[m,] = NA
    current_gse_featureData[m,] = NA}
}
current_gse_assayData = na.omit(current_gse_assayData)
current_gse_featureData = na.omit(current_gse_featureData)
current_gse_assayData = data.frame(current_gse_assayData)
#assigning ENTREZIDs to current_gse_assayData
rownames(current_gse_assayData) = unlist(current_gse_featureData["ENTREZID"])
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
color_phd = paste(current_gse_phenoData[c(1:6,19:24),2],current_gse_phenoData[c(1:6,19:24),3],sep="/")
cluster_plot = cluster_plot + geom_point(aes(color = color_phd))
cluster_plot

#limma
top_genes_limma = ""
column_subset = c(1:3,4:6)
#create assayData table
gse_assayData_for_limma = current_gse_assayData[,column_subset]
#create des_matrix
desmatrix = cbind(c(rep(1,6)), c(c(rep(1,3)), c(rep(0,3))))
colnames(desmatrix) = c("intercept", "phenotype")
rownames(desmatrix) = colnames(gse_assayData_for_limma)
#create cont_matrix
cont_matrix = matrix(c(0,1), nrow=2, ncol = 1)

limma_maker(gse_assayData_for_limma, desmatrix, cont_matrix)
View(top_genes_limma)  
top_genes_control = top_genes_limma


#correlation
top_genes_all <- top_genes_all[order(row.names(top_genes_all)),]
top_genes_5month <- top_genes_5month[order(row.names(top_genes_5month)),]
top_genes_24month <- top_genes_24month[order(row.names(top_genes_24month)),]
top_genes_control <- top_genes_control[order(row.names(top_genes_control)),]

un5_24 = union(rownames(top_genes_5month[top_genes_5month$adj.P.Val<0.05,]), rownames(top_genes_24month[top_genes_24month$adj.P.Val<0.05,]))
unall_5 = union(rownames(top_genes_all[top_genes_all$adj.P.Val<0.05,]), rownames(top_genes_5month[top_genes_5month$adj.P.Val<0.05,]))
unall_24 = union(rownames(top_genes_all[top_genes_all$adj.P.Val<0.05,]), rownames(top_genes_24month[top_genes_24month$adj.P.Val<0.05,]))
top_genes_all_for_5 = subset(top_genes_all, rownames(top_genes_all) %in% unall_5)
top_genes_5_for_all = subset(top_genes_5month, rownames(top_genes_5month) %in% unall_5)
top_genes_all_for_24 = subset(top_genes_all, rownames(top_genes_all) %in% unall_24)
top_genes_24_for_all = subset(top_genes_24month, rownames(top_genes_24month) %in% unall_24)
top_genes_24_for_5 = subset(top_genes_24month, rownames(top_genes_24month) %in% un5_24)
top_genes_5_for_24 = subset(top_genes_5month, rownames(top_genes_5month) %in% un5_24)

#DRAFT
#create vennDiagramm for interception of genes

top_genes_all_suff_pval = top_genes_all[top_genes_all$adj.P.Val<0.05,]
top_genes_5month_suff_pval = top_genes_5month[top_genes_5month$adj.P.Val<0.05,]
top_genes_24month_suff_pval = top_genes_24month[top_genes_24month$adj.P.Val<0.05,]
top_genes_control_suff_pval = top_genes_control[top_genes_control$adj.P.Val<0.05,]

top_genes_all_plus = top_genes_all_suff_pval[top_genes_all_suff_pval$logFC>0,]
top_genes_all_minus = top_genes_all_suff_pval[top_genes_all_suff_pval$logFC<0,]
top_genes_5month_plus = top_genes_5month_suff_pval[top_genes_5month_suff_pval$logFC>0,]
top_genes_5month_minus = top_genes_5month_suff_pval[top_genes_5month_suff_pval$logFC<0,]
top_genes_24month_plus = top_genes_24month_suff_pval[top_genes_24month_suff_pval$logFC>0,]
top_genes_24month_minus = top_genes_24month_suff_pval[top_genes_24month_suff_pval$logFC<0,]
top_genes_control_plus = top_genes_control_suff_pval[top_genes_control_suff_pval$logFC>0,]
top_genes_control_minus = top_genes_control_suff_pval[top_genes_control_suff_pval$logFC<0,]

#diffexpressed genes with sufficient p-val
length(rownames(top_genes_all[top_genes_all$adj.P.Val<0.05,]))
length(rownames(top_genes_5month[top_genes_5month$adj.P.Val<0.05,]))
length(rownames(top_genes_24month[top_genes_24month$adj.P.Val<0.05,]))


grid.newpage()
draw.triple.venn(area1 = length(rownames(top_genes_all_plus)), 
                 area2 = length(rownames(top_genes_5month_plus)),
                 area3 = length(rownames(top_genes_24month_plus)),
                 n12 = length(intersect(rownames(top_genes_all_plus), rownames(top_genes_5month_plus))),
                 n23 = length(intersect(rownames(top_genes_5month_plus), rownames(top_genes_24month_plus))),
                 n13 = length(intersect(rownames(top_genes_all_plus), rownames(top_genes_24month_plus))), 
                 n123 = length(intersect(intersect(rownames(top_genes_all_plus), rownames(top_genes_5month_plus)), rownames(top_genes_24month_plus))), 
                 category = c("All", "5 month", "24 month"), lty = "blank", 
                 fill = c("skyblue", "pink", "magenta"))

grid.newpage()
draw.triple.venn(area1 = length(rownames(top_genes_all_minus)), 
                 area2 = length(rownames(top_genes_5month_minus)),
                 area3 = length(rownames(top_genes_24month_minus)),
                 n12 = length(intersect(rownames(top_genes_all_minus), rownames(top_genes_5month_minus))),
                 n23 = length(intersect(rownames(top_genes_5month_minus), rownames(top_genes_24month_minus))),
                 n13 = length(intersect(rownames(top_genes_all_minus), rownames(top_genes_24month_minus))), 
                 n123 = length(intersect(intersect(rownames(top_genes_all_minus), rownames(top_genes_5month_minus)), rownames(top_genes_24month_minus))), 
                 category = c("All", "5 month", "24 month"), lty = "blank", 
                 fill = c("skyblue", "pink", "magenta"))

grid.newpage()
draw.triple.venn(area1 = length(rownames(top_genes_control_minus)), 
                 area2 = length(rownames(top_genes_5month_plus)),
                 area3 = length(rownames(top_genes_24month_plus)),
                 n12 = length(intersect(rownames(top_genes_control_minus), rownames(top_genes_5month_plus))),
                 n23 = length(intersect(rownames(top_genes_5month_plus), rownames(top_genes_24month_plus))),
                 n13 = length(intersect(rownames(top_genes_control_minus), rownames(top_genes_24month_plus))), 
                 n123 = length(intersect(intersect(rownames(top_genes_control_minus), rownames(top_genes_5month_plus)), rownames(top_genes_24month_plus))), 
                 category = c("control", "5 month", "24 month"), lty = "blank", 
                 fill = c("skyblue", "pink", "magenta"))

grid.newpage()
draw.triple.venn(area1 = length(rownames(top_genes_control_plus)), 
                 area2 = length(rownames(top_genes_5month_minus)),
                 area3 = length(rownames(top_genes_24month_minus)),
                 n12 = length(intersect(rownames(top_genes_control_plus), rownames(top_genes_5month_minus))),
                 n23 = length(intersect(rownames(top_genes_5month_minus), rownames(top_genes_24month_minus))),
                 n13 = length(intersect(rownames(top_genes_control_plus), rownames(top_genes_24month_minus))), 
                 n123 = length(intersect(intersect(rownames(top_genes_control_plus), rownames(top_genes_5month_minus)), rownames(top_genes_24month_minus))), 
                 category = c("control", "5 month", "24 month"), lty = "blank", 
                 fill = c("skyblue", "pink", "magenta"))
