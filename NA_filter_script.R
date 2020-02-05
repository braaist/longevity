#NA filtering for datasets
union_rownames = c()
genes_in_datasets = data.frame()
for (i in gse_list_main){
  union_rownames= union(rownames(i), union_genes)
}
logFCmatrix = matrix(nrow = length(union_rownames), ncol = length(gse_list_main))
rownames(logFCmatrix) = union_rownames
colnames(logFCmatrix) = names(gse_list_main)
for (i in names(gse_list_main)){
  for (j in rownames(gse_list_main[[i]])){
    logFCmatrix[j,i] = gse_list_main[[i]][j,1]
  }
}

logFCmatrix$NACount = rowSums(is.na(logFCmatrix))
ggplot(logFCmatrix, aes(x = NACount)) + geom_density()