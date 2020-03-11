library(reshape2)
library(psych)
library(BlandAltmanLeh)
library(igraph)
library(tidyverse)
library(reshape2)


#create logFCmatrix main
rownames_vec = c()
colnames_vec = c()
for (i in gse_list_main){
  rownames_vec = c(rownames_vec, rownames(i))
}
rownames_vec = unique(rownames_vec)
colnames_vec = names(gse_list_main)
logFCmatrix = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(logFCmatrix) = rownames_vec
colnames(logFCmatrix) = colnames_vec
for (i in rownames_vec){
  for (name in names(gse_list_main)){
    current_ds = data.frame(gse_list_main[name])
    if (i %in% rownames(current_ds)){
      FCname = paste(name, ".logFC", sep = "")
      logFCmatrix[i,name] = current_ds[i, 1]
    }
  }
}

#create SE matrix
rownames_vec = c()
colnames_vec = c()
for (i in gse_list_main){
  rownames_vec = c(rownames_vec, rownames(i))
}
rownames_vec = unique(rownames_vec)
colnames_vec = names(gse_list_main)
SEmatrix = matrix(nrow = length(rownames_vec), ncol = length(colnames_vec))
rownames(SEmatrix) = rownames_vec
colnames(SEmatrix) = colnames_vec
for (i in rownames_vec){
  for (name in names(gse_list_main)){
    current_ds = data.frame(gse_list_main[name])
    if (i %in% rownames(current_ds)){
      FCname = paste(name, ".logFC", sep = "")
      SEmatrix[i,name] = current_ds[i, 1]
    }
  }
}



#cormatrixsign
cormatrixsign = list()
cormatrixsign[["100"]] = data.frame()
cormatrixsign[["200"]] = data.frame()
cormatrixsign[["300"]] = data.frame()
cormatrixsign[["400"]] = data.frame()
cormatrixsign[["500"]] = data.frame()
cormatrixsign[["750"]] = data.frame()
cormatrixsign[["1000"]] = data.frame()
cormatrixsign[["1500"]] = data.frame()
cormatrixsign[["2000"]] = data.frame()
for (thres in names(cormatrixsign)){
  for (i in 1:length(gse_list_main)){
    for (j in i:length(gse_list_main)){
      topA = gse_list_main[[i]] %>% rownames_to_column(var = "row.names")
      topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topA = topA %>% column_to_rownames(var = "row.names")
      topB = gse_list_main[[j]] %>% rownames_to_column(var = "row.names")
      topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topB = topB %>% column_to_rownames(var = "row.names")
      totalrownames = union(rownames(topA), rownames(topB))
      cormatrixsign[[thres]][names(gse_list_main)[i], names(gse_list_main)[j]] = cor(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = "spearman", use = "complete.obs")
      cormatrixsign[[thres]][names(gse_list_main)[j], names(gse_list_main)[i]] = cor(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = "spearman", use = "complete.obs")
      #    mergedmatrix = gse_list_main[[i]]["logFC"]
      #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
      #    mergedmatrix = merge(mergedmatrix, gse_list_main[[j]]["logFC"], by=0, all=TRUE)
      #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
      #    cormatrixdenoised[names(gse_list_main)[i], names(gse_list_main)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
    }
  }
}

cormatrixsign[["all"]] = as.data.frame(cor(logFCmatrix, method = "spearman", use = "pairwise.complete.obs"))


corpvalsign = list()
corpvalsign[["100"]] = data.frame()
corpvalsign[["200"]] = data.frame()
corpvalsign[["300"]] = data.frame()
corpvalsign[["400"]] = data.frame()
corpvalsign[["500"]] = data.frame()
corpvalsign[["750"]] = data.frame()
corpvalsign[["1000"]] = data.frame()
corpvalsign[["1500"]] = data.frame()
corpvalsign[["2000"]] = data.frame()
for (thres in names(corpvalsign)){
  for (i in 1:length(gse_list_main)){
    for (j in i:length(gse_list_main)){
      topA = gse_list_main[[i]] %>% rownames_to_column(var = "row.names")
      topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topA = topA %>% column_to_rownames(var = "row.names")
      topB = gse_list_main[[j]] %>% rownames_to_column(var = "row.names")
      topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topB = topB %>% column_to_rownames(var = "row.names")
      totalrownames = union(rownames(topA), rownames(topB))
      corpvalsign[[thres]][names(gse_list_main)[i], names(gse_list_main)[j]] = as.numeric(cor.test(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = "spearman")$p.value)
      corpvalsign[[thres]][names(gse_list_main)[j], names(gse_list_main)[i]] = as.numeric(cor.test(gse_list_main[[i]][totalrownames,]$logFC, gse_list_main[[j]][totalrownames,]$logFC, method = "spearman")$p.value)
      #    mergedmatrix = gse_list_main[[i]]["logFC"]
      #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
      #    mergedmatrix = merge(mergedmatrix, gse_list_main[[j]]["logFC"], by=0, all=TRUE)
      #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
      #    cormatrixdenoised[names(gse_list_main)[i], names(gse_list_main)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
    }
  }
}

corpvalsign[["all"]] = data.frame()
for (colnum1 in 1:length(colnames(logFCmatrix))){
  for (colnum2 in colnum1:length(colnames(logFCmatrix))){
    corpvalsign[["all"]][colnames(logFCmatrix)[colnum1], colnames(logFCmatrix)[colnum2]] = as.numeric(cor.test(logFCmatrix[, colnum1], logFCmatrix[, colnum2], method = "spearman")$p.value)
    corpvalsign[["all"]][colnames(logFCmatrix)[colnum2], colnames(logFCmatrix)[colnum1]] = as.numeric(cor.test(logFCmatrix[, colnum1], logFCmatrix[, colnum2], method = "spearman")$p.value)
  }
}

coradjpvalsign = list()
coradjpvalsign[["100"]] = data.frame()
coradjpvalsign[["200"]] = data.frame()
coradjpvalsign[["300"]] = data.frame()
coradjpvalsign[["400"]] = data.frame()
coradjpvalsign[["500"]] = data.frame()
coradjpvalsign[["750"]] = data.frame()
coradjpvalsign[["1000"]] = data.frame()
coradjpvalsign[["1500"]] = data.frame()
coradjpvalsign[["2000"]] = data.frame()
coradjpvalsign[["all"]] = data.frame()
for (thres in names(corpvalsign)){
  vec = as.vector(corpvalsign[[thres]][upper.tri(corpvalsign[[thres]], diag = F)])
  vec = p.adjust(vec, method = "BH")
  tempmatrix = matrix(0, length(gse_list_main), length(gse_list_main))
  tempmatrix[lower.tri(tempmatrix, diag = F)] = vec
  tempmatrix = t(tempmatrix)
  tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
  coradjpvalsign[[thres]] = tempmatrix
  colnames(coradjpvalsign[[thres]]) = colnames(corpvalsign[[thres]])
  rownames(coradjpvalsign[[thres]]) = rownames(corpvalsign[[thres]])
}

corofcor = matrix(nrow = length(names(cormatrixsign)), ncol = length(names(cormatrixsign)))
colnames(corofcor) = names(cormatrixsign)
rownames(corofcor) = names(cormatrixsign)
for (thres1 in names(cormatrixsign)){
  # vectorize the upper triangles of matrices (diagonals excluded)
  vec1 = as.vector(as.matrix(cormatrixsign[[thres1]])[upper.tri(as.matrix(cormatrixsign[[thres1]]), diag = F)])
  for (thres2 in names(cormatrixsign)){
    vec2 = as.vector(as.matrix(cormatrixsign[[thres2]])[upper.tri(as.matrix(cormatrixsign[[thres2]]), diag = F)])
    corofcor[thres1, thres2] = cor(vec1, vec2, method = "spearman")
  }
}
# View the result
View(corofcor)

cortestsign = list()
cortestsign[["100"]] = data.frame()
cortestsign[["200"]] = data.frame()
cortestsign[["300"]] = data.frame()
cortestsign[["400"]] = data.frame()
cortestsign[["500"]] = data.frame()
cortestsign[["750"]] = data.frame()
cortestsign[["1000"]] = data.frame()
cortestsign[["1500"]] = data.frame()
cortestsign[["2000"]] = data.frame()
cortestsign[["all"]] = data.frame()
for (thres in names(cormatrixsign)){
  for (colname in colnames(cormatrixsign[[thres]])){
    for (rowname in rownames(cormatrixsign[[thres]])){
      if (coradjpvalsign[[thres]][rowname, colname] < 0.05){
        if (cormatrixsign[[thres]][rowname, colname] > 0){
          cortestsign[[thres]][rowname, colname] = 1
        } else {
          cortestsign[[thres]][rowname, colname] = -1
        } 
      } else{
        cortestsign[[thres]][rowname, colname] = 0
      }
    }
  }
}


for (thres in names(cormatrixsign)){
  print(shapiro.test(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F)]))
  qqnorm(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F)])
  qqline(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F)])
}


for (thres1 in names(cormatrixsign)){
  tempmatrix = matrix(0, length(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]), 2)
  tempmatrix[,1] = cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]
  for (thres2 in names(cormatrixsign)){
    if (thres1 == thres2) {
      break
    }
    tempmatrix[,2] = cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]
    tempmatrix = as.data.frame(tempmatrix)
    colnames(tempmatrix) = c(paste0("threshold_", thres1), paste0("threshold_", thres2))
    a = ggplot(tempmatrix, aes_string(x = colnames(tempmatrix)[1], y = colnames(tempmatrix)[2])) + geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_abline(intercept = 0)
    print(a)
  }
}


# Bland-Altman plots for absolute value of correlation:
utestpval = matrix(nrow = length(names(cormatrixsign)), ncol = length(names(cormatrixsign)))
colnames(utestpval) = names(cormatrixsign)
rownames(utestpval) = names(cormatrixsign)
for (thres1 in names(cormatrixsign)){
  for (thres2 in names(cormatrixsign)){
    if (thres1 == thres2){
      break
    }
    meandif = round(mean(abs(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]) - abs(cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)])), 5)
    a = as.numeric(wilcox.test(abs(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]), abs(cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]))$p.value)
    utestpval[thres1, thres2] = a
    utestpval[thres2, thres1] = a
    bland.altman.plot(abs(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]), abs(cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]), main=paste0("Bland Altman Plot of ", thres1, " against ", thres2, " (Mann-Wintey U test p-value: ", a, "; mean difference: ", meandif, ")"), xlab="Means", ylab=paste0("Differences (", thres1, " - ", thres2, ")"))
  }
}

# calculate mean absolute correlation values for all thresholds and plot it against correlations with the "no threshold" correlation set
meanabscor = vector()
i = 1
for (thres in names(cormatrixsign)){
  meanabscor[i] = mean(abs(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F) & cortestsign[[thres]] != 0]))
  i = i + 1
}
tableforggplot = as.data.frame(corofcor[,"all"])
tableforggplot$meanabscor = meanabscor
tableforggplot$index = rownames(tableforggplot)
colnames(tableforggplot) = c("corwithall", "meanabscor", "index")
#tableforggplot$index = factor(tableforggplot$index, levels = tableforggplot$index)
tableforggplot$index = as.numeric(as.character(tableforggplot$index))
meltedtableforggplot = melt(tableforggplot, id = "index")
ggplot(meltedtableforggplot, aes(x = index, y = value, colour = variable)) + geom_point()

portiontable = matrix(nrow = length(names(cortestsign)), ncol = length(names(cortestsign)))
colnames(portiontable) = names(cortestsign)
rownames(portiontable) = names(cortestsign)
for (thres1 in names(cortestsign)){
  for (thres2 in names(cortestsign)){
    consistentones = sum(cortestsign[[thres1]][upper.tri(cortestsign[[thres1]], diag = F)] == cortestsign[[thres2]][upper.tri(cortestsign[[thres2]], diag = F)])
    total = length(cortestsign[[thres1]][upper.tri(cortestsign[[thres1]], diag = F)])
    portiontable[thres1, thres2] = consistentones / total
    portiontable[thres2, thres1] = consistentones / total
  }
}
# check out the result:
View(portiontable)
portiontable = data.frame(portiontable)
ggplot(portiontable, aes(x = index, y = value, colour = variable)) + geom_point()

View(meltedtableforggplot)


# distribution of correlation coefficients:
chosencormatrix = cormatrixsign[["500"]]
chosenpvalmatrix = coradjpvalsign[["500"]]
chosentestmatrix = cortestsign[["500"]]
wilcox.test(chosencormatrix[upper.tri(chosencormatrix, diag = F)])

# cor distribution:
matrixforggplot = as.data.frame(chosencormatrix[upper.tri(chosencormatrix, diag = F)])
colnames(matrixforggplot) = c("value")
ggplot(matrixforggplot, aes(x = value)) + geom_density()+geom_vline(xintercept = 0) + geom_vline(xintercept = median(as.numeric(matrixforggplot$value)),  color = "red", linetype = "dashed")



cormat = cormatrixsign[[5]]
corpvalmat = coradjpvalsign[[5]]

cormatrix = apply(cormat, 2, rev)
upper_tri <- get_upper_tri(cormatrix)
# Melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat <- melt(cormatrix, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggheatmap


pops = cormat>0.1 & corpvalmat<0.05
pops = data.frame(pops)
sum_vec_plus = pops %>% summarise_all(funs(sum))
sum_vec_plus = data.frame(t(sum_vec_plus))
sum_vec_plus$rn =  rownames(sum_vec_plus)
sum_vec_plus$name = "plus"
colnames(sum_vec_plus) = c("number", "name", "status")

p<-ggplot(data=sum_vec, aes(x=rownames(sum_vec_plus), y=c(sum_vec_plus$t.sum_vec_plus., sum_vec_minus$t.sum_vec_minus.))) +
  geom_bar(stat="identity", fill = "red")
p

flops = cormat< -0.1 & corpvalmat<0.05
flops = data.frame(flops)
sum_vec_minus = flops %>% summarise_all(funs(sum))
sum_vec_minus = data.frame(t(sum_vec_minus))
sum_vec_minus$rn = rownames(sum_vec_minus)
sum_vec_minus$name = "minus"
colnames(sum_vec_minus) = c("number", "name", "status")


sum_vec_all = bind_rows(sum_vec_plus, sum_vec_minus)


test = cormat< -0.1

p <- sum_vec_all %>% ggplot(aes(x = number, fill = status)) + geom_bar(stat="identity", aes(y = rownames(sum_vec_all)))
p

new.df = c()
new.df<-melt(cortestsign[["500"]]) 
new.df<-new.df[complete.cases(new.df),]
new.df$Score<-factor(new.df$value, levels = c("1", "-1", "0"))
new.df = new.df %>% filter(Score != 0)
ggplot(new.df, aes(x=variable, fill=Score)) + geom_bar()

ggplot(sum_vec_all, aes(fill=status, y=number, x=name)) + geom_bar(position="stack", stat="identity") + 
  scale_fill_manual("legend", values = c("minus" = "blue", "plus" = "red"))


