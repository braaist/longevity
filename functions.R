entrez_converter <- function(ACCNUM){
  #function for convertation ACCNUM to ENTREZID
  if(is.null(xx[[as.character(ACCNUM)]])){
    return(NA)
  }else{
    return(xx[[as.character(ACCNUM)]])}
}

limma_maker <- function(assayData, des_matrix, cont_matrix){
  #function for getting data about differential expression and creating volcano plots
  fit = lmFit(assayData, des_matrix)
  fit_contrast = eBayes(contrasts.fit(fit, cont_matrix))
  volcanoplot(fit_contrast)
  top_genes = topTable(fit_contrast, number = Inf, adjust = "BH")
  top_genes_limma <<- top_genes
  return(top_genes)
}
