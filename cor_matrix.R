#SCRIPT for calculating correlation matrix and building heatmap from list with GSE

library(ggplot2)
library(reshape2)


cormat = data.frame()
for (i in names(gse_list_main)){
  for (j in names(gse_list_main)){
    pops = (union(rownames(gse_list_main[[i]])[1:150], rownames(gse_list_main[[j]])[1:150]))
    temp_m = data.frame()
    for (m in pops){
      temp_m[m,1] = gse_list_main[[i]][m,1]
      temp_m[m,2] = gse_list_main[[j]][m,1]
    }
    print(cor(x = temp_m[,1], y = temp_m[,2], use = "complete.obs", method = "spearman"))
    cormat[i,j] = cor(x = temp_m[,1], y = temp_m[,2], use = "complete.obs", method = "spearman")
  }
}


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Reorder the correlation matrix
#cormatrix <- reorder_cormat(cormat)
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


