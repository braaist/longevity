#get GSE file create phenodata table
library(GEOquery)
current_gse <- getGEO("gse55272")
#create phenodata table
current_gse_phenodata = data.frame(pData(phenoData(current_gse[[1]]))[,c(44,41,40,42,43,25)])
current_gse_phenodata = subset(gse55272_phd, tissue.ch1=="Liver")
current_gse_phenodata = current_gse_phenodata %>% mutate("diet" = "-", "drug" = "-", "dose" = "-","age_intervention" = "-")
#create assayData table
current_gse_assayData = assayData((current_gse[[1]]))[["exprs"]]
