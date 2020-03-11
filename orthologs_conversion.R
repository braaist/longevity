library(biomaRt)
library(tidyverse)
library(tools)
#SASHA'S script for orthologs conversion

#gse120081 - human
#gse28688
#gse81891
#gse89455

#convert Rda to Rds
for (i in list.files("Users/braaist/Documents/longevity/longevity/gse127927/")){
  if (file_ext(i) == 'Rda'){
    a = paste(file_path_sans_ext(i), ".Rds", sep="")
    load(paste("Users/braaist/Documents/longevity/longevity/gse127927/", as.character(i), sep = ""))
    saveRDS(top_genes_limma, file = paste("Users/braaist/Documents/longevity/longevity/gse127927/", as.character(a), sep = "") )
  }
}



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
human_mouse_entrez_map = na.omit(Human_to_mouse_orthologs1)
human_mouse_entrez_map = as.data.frame(human_mouse_entrez_map)
human_mouse_entrez_map = human_mouse_entrez_map %>% mutate_all(as.character)

#gse28688
flops = readRDS("gse28688/top_genes_HFF_gse28688.Rds")
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

saveRDS(flops, file="gse28688/ort_top_genes_HFF_gse28688.Rds")

#gse81891
flops = readRDS("gse81891/top_genes_OSKM_gse81891.Rds")
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

saveRDS(flops, file="gse81891/ort_top_genes_OSKM_gse81891.Rds")



#gse89455
for (name in list.files("Users/braaist/Documents/longevity/longevity/gse89455/")){
  if (file_ext(name) == "Rds"){
    print(name)

      flops = readRDS(paste("Users/braaist/Documents/longevity/longevity/gse89455/", name, sep=""))
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
        if ((length(which(flops[i,'orthologs'] == flops[,'orthologs'])) != 1) && (list(which(flops[i,'orthologs'] == flops[,'orthologs'])) %notin% same_ort_list)){
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
      flops_out = paste("ort_", name, sep = "")
      print(flops_out)
      file = paste("Users/braaist/Documents/longevity/longevity/gse89455/", flops_out, sep="")
      saveRDS(flops, file)
    }
  }
}


#final human signature conversion
flops = human_signature
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

human_ort_signature = flops




