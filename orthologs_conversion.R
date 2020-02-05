library(biomaRt)
library(tidyverse)
#SASHA'S script for orthologs conversion

#gse120081 - human
#gse28688
#gse81891
#gse89455

ensembl = useEnsembl("ensembl", host = "uswest.ensembl.org")

datasets <- listDatasets(ensembl)
human_dataset = useDataset("hsapiens_gene_ensembl", mart=ensembl)
mouse_dataset = useDataset("mmusculus_gene_ensembl", mart=ensembl)

human = useEnsembl("ensembl", "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
mouse = useEnsembl("ensembl","mmusculus_gene_ensembl", host = "uswest.ensembl.org")

Human_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                   mart=human_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Human_to_mouse_orthologs) <- c("Human_Entrez","Mouse_Entrez")

Human_to_mouse_orthologs1 = Human_to_mouse_orthologs %>% group_by(Human_Entrez) %>% filter(n() == 1)
Human_to_mouse_orthologs1 = na.omit(Human_to_mouse_orthologs1)
Mouse_to_human_orthologs1 = subset(Mouse_to_human_orthologs, Mouse_Entrez %in% Human_to_mouse_orthologs1$Mouse_Entrez)
human_mouse_entrez_map = na.omit(Human_to_mouse_orthologs1)
human_mouse_entrez_map = as.data.frame(human_mouse_entrez_map)
human_mouse_entrez_map = human_mouse_entrez_map %>% mutate_all(as.character)

#gse28688
flops = readRDS("all_FClists2/top_genes_HFF_gse28688.Rds")
View(flops)
flops$orthologs = NA
flops$original = rownames(flops)
rownames(flops) = 1:nrow(flops)
for (i in rownames(flops)){
  if (flops[i, "original"] %in% human_mouse_entrez_map$Human_Entrez){
    flops[i, "orthologs"] = human_mouse_entrez_map[human_mouse_entrez_map$Human_Entrez == flops[i,"original"],"Mouse_Entrez"]
  }
}

NAfilter = !is.na(flops$orthologs)
flops = flops[NAfilter,]
rownames(flops) = flops$orthologs


same_ort_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(flops)){
  print(i)
  if ((length(which(flops[i,'orthologs'] == flops[,'orthologs'])) != 1) && (list(which(flops[i,'orthologs'] == flops[,'orthologs'])) %notin% same_entrez_list)){
    same_ort_list = append(same_ort_list, list(which(flops[i,'orthologs'] == flops[,'orthologs'])))
  }
}

for (i in same_ort_list){
  i = i[-1]
  for (m in i){
    flops[m,] = NA
  }
}

flops = na.omit(flops)
rownames(flops) = flops$orthologs
flops$original = NULL
flops$orthologs = NULL

saveRDS(flops, file="all_FClists2/top_genes_HFF_gse28688_ort.Rds")

#gse81891
flops = readRDS("all_FClists2/top_genes_OS(v)KM_gse81891.Rds")
flops = flops[1:1000,]
flops$orthologs = NA
flops$original = rownames(flops)
rownames(flops) = 1:nrow(flops)
for (i in rownames(flops)){
  if (flops[i, "original"] %in% human_mouse_entrez_map$Human_Entrez){
    flops[i, "orthologs"] = human_mouse_entrez_map[human_mouse_entrez_map$Human_Entrez == flops[i,"original"],"Mouse_Entrez"]
  }
}

NAfilter = !is.na(flops$orthologs)
flops = flops[NAfilter,]


same_ort_list = list()
'%notin%' = Negate('%in%')
for (i in rownames(flops)){
  print(i)
  if ((length(which(flops[i,'orthologs'] == flops[,'orthologs'])) != 1) && (list(which(flops[i,'orthologs'] == flops[,'orthologs'])) %notin% same_entrez_list)){
    same_ort_list = append(same_ort_list, list(which(flops[i,'orthologs'] == flops[,'orthologs'])))
  }
}

for (i in same_ort_list){
  i = i[-1]
  for (m in i){
    flops[m,] = NA
  }
}

flops = na.omit(flops)
rownames(flops) = flops$orthologs
flops$original = NULL
flops$orthologs = NULL

saveRDS(flops, file="all_FClists2/top_genes_OS(v)KM_gse81891_ort.Rds")



#gse89455
for (name in list.files("all_FClists2/")){
  a = strsplit(name, "_")
  for (m in a){
    if(m[4] == "gse89455.Rds"){
      
      flops = readRDS(paste("all_FClists2/", name, sep=""))
      flops = flops[1:1000,]
      flops$orthologs = NA
      flops$original = rownames(flops)
      rownames(flops) = 1:nrow(flops)
      for (i in rownames(flops)){
        if (flops[i, "original"] %in% human_mouse_entrez_map$Human_Entrez){
          flops[i, "orthologs"] = human_mouse_entrez_map[human_mouse_entrez_map$Human_Entrez == flops[i,"original"],"Mouse_Entrez"]
        }
      }
      
      NAfilter = !is.na(flops$orthologs)
      flops = flops[NAfilter,]
      
      
      same_ort_list = list()
      '%notin%' = Negate('%in%')
      for (i in rownames(flops)){
        print(i)
        if ((length(which(flops[i,'orthologs'] == flops[,'orthologs'])) != 1) && (list(which(flops[i,'orthologs'] == flops[,'orthologs'])) %notin% same_entrez_list)){
          same_ort_list = append(same_ort_list, list(which(flops[i,'orthologs'] == flops[,'orthologs'])))
        }
      }
      
      for (i in same_ort_list){
        i = i[-1]
        for (m in i){
          flops[m,] = NA
        }
      }
      
      flops = na.omit(flops)
      rownames(flops) = flops$orthologs
      flops$original = NULL
      flops$orthologs = NULL
      flops_out = paste("ort_cut_", name, sep = "")
      print(flops_out)
      saveRDS(flops, file=paste("all_FClists2/", flops_out, sep="", ))
    }
  }
}






