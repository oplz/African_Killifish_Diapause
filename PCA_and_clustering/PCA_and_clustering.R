## PCA and cluster plotting
##

library(FactoMineR) 
library(pheatmap)

## working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory


##### ===========================================================================================
##### data.frame making 
##### ===========================================================================================

## to make a list of all libraries in the data dicretory 
datadir <- dir("read_mapping/")

## to create a matrix to include data from all libraries 
n <- length(datadir)

lengthmeasure<-read.table(paste0("read_mapping/01/abundance.tsv"), header=T)
m <- length(lengthmeasure[,1])

tpm.matrix <- matrix(numeric(length = n * m),
                     ncol = n,
                     nrow = m)
colnames(tpm.matrix) <- datadir

for( i in 1:n ){
    data.i <- read.table(paste0("read_mapping/",
                                datadir[i],
                                "/abundance.tsv"),
                         colClasses = c("character", "NULL","NULL", "numeric", "numeric"),
                         sep = "\t",
                         header = TRUE)
    tpm.matrix[, i] <- data.i$tpm
    tsvref<-data.i
    rm(data.i)
    cat(i, "\n")
}

rownames(tpm.matrix)<-tsvref$target_id  

save(tpm.matrix,
     file = "PCA_and_Clustering/output/post_combine_kallisto_results_matrices.RData")


## to rotate matrix: PCA and many other functions require individual in rows and genes in columns
rtpm <- t(apply(tpm.matrix, 2, rev))

## to assign library names
rownames(rtpm)<-c("02-Pre", "03-Pre", "01-Pre", 
                  "04-Non", "05-Non", "06-Non", 
                  "08-Dia1m", "10-Dia1m", "09-Dia1m", 
                  "11-Dia6d", "12-Dia6d", "13-Dia6d", 
                  "14-Dia3d", "15-Dia3d", "16-Dia3d")

rm(list=setdiff(ls(), "rtpm"))


##### ===========================================================================================
##### PCA plotting
##### ===========================================================================================

pdf("PCA_and_Clustering/output/post_combine_PCA.pdf",height=9,width=9)

pca1=PCA(log2(rtpm+1), graph=F)

## to reverse PC3 axis to make it consistent with PC2 (from pre-diapause to non-diapause)
pca1$ind$coord[,3]<-pca1$ind$coord[,3]*-1

datacolor<-as.vector(c("#FFCC99", "#FFCC99", "#FFCC99", 
                       "#F48000", "#F48000", "#F48000",
                       "#0070C0", "#0070C0", "#0070C0",
                       "#5B9BD5", "#5B9BD5", "#5B9BD5",
                       "#9BC2E6", "#9BC2E6", "#9BC2E6"))

plot(pca1$eig[,2], type="p", col="blue", pch=16,
     main="Plot of Principle Components", 
     xlab="Principle Components", 
     ylab="Percentage of Variance")

plot.PCA(pca1, 
         choix="ind",
         axes=c(1,2),
         habillage="ind", 
         label="none", 
         cex=4,
         title="Plot of Individuals (PC1 vs PC2)",
         palette=palette(datacolor)
)

plot.PCA(pca1, 
         choix="ind",
         axes=c(2,3),
         habillage="ind", 
         label="none", 
         cex=4,
         title="Plot of Individuals (PC2 vs PC3)",
         palette=palette(datacolor)
)

plot.PCA(pca1, 
         choix="ind",
         axes=c(1,3),
         habillage="ind", 
         label="none", 
         cex=4,
         title="Plot of Individuals (PC1 vs PC3)",
         palette=palette(datacolor)
)



## to extract contribution of all variables to PC Dim 1 - 5
dimtank<-as.matrix(pca1$var$contrib)
write.csv(dimtank, 
          file = "PCA_and_Clustering/output/post_combine_PC1to5_gene_contribution.csv")

for(i in 1:5) 
{ 
  handler1<-as.matrix(dimtank[,i])
  rownames(handler1)<-rownames(dimtank)
  handler2<-sort(handler1[,1], decreasing=T)
  nam <- paste("contributor", i, sep = "")
  assign(nam, as.matrix(handler2, row.names=rownames(handler2)))
  rm(handler1, handler2, nam)
}

save(contributor1,contributor2, contributor3, contributor4, contributor5, 
     file = "PCA_and_Clustering/output/post_combine_genes_of_PC.RData")

dev.off()

rm(list=setdiff(ls(), "rtpm"))


##### ===========================================================================================
##### Hierachical Clustering % Heatmap
##### ===========================================================================================

pdf("PCA_and_Clustering/output/post_combined_clustering.pdf",height=9,width=9)

## clustering
datamat <- as.matrix(rtpm)
datalog <- log2(datamat+1)
Eucildean_distance <- dist(scale(datalog, center=T, scale=T))   # compute the dissimilarity matrix of the standardised data using Eucildean distance
clust_results <- hclust(Eucildean_distance, method="ward.D2")   # calculate a hieracharchical clustering of data
plot(clust_results, sub=c(""))

## heatmap
color_panel<-colorRampPalette(c("red","lightyellow","darkblue"))
orddis<-as.matrix(Eucildean_distance)
pheatmap(orddis, symm=T, main="the original dissimilarity matrix", color=color_panel(100))

dev.off()

rm(list=ls())

