#get GSE file create phenodata table
library(GEOquery)
library(dplyr)
library(org.Mm.eg.db)
current_gse = getGEO("gse55272")
#create phenodata table
current_gse_phenodata = data.frame(pData(phenoData(current_gse[[1]]))[,c(44,41,40,42,43,25)])
current_gse_phenodata = subset(gse55272_phd, tissue.ch1=="Liver")
current_gse_phenodata = current_gse_phenodata %>% mutate("diet" = "-", "drug" = "-", "dose" = "-","age_intervention" = "-")
#create assayData table
current_gse_assayData = assayData((current_gse[[1]]))[["exprs"]]
#Making lognorm plots
lognorm_current_gse_assayData = scale(log(current_gse_assayData))
den = apply(lognorm_current_gse_assayData, 2, density)
plot(NA, xlim=range(sapply(den, "[", "x")), ylim=range(sapply(den, "[", "y")))
mapply(lines, den, col=1:length(den))