## This script is for Nfur embryogenesis and diapause RNA-Seq library
##

library(FactoMineR) 
library(pheatmap)

setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

##### ===========================================================================================
##### data.frame making 
##### ===========================================================================================


## to make a list of all libraries in the data dicretory 
datadir <- dir("se_read_mapping/")

## to create a matrix to include data from all libraries 
n <- length(datadir)

lengthmeasure<-read.table(paste0("se_read_mapping/0D239/abundance.tsv"), header=T)
m <- length(lengthmeasure[,1])

tpm.matrix <- matrix(numeric(length = n * m),
                     ncol = n,
                     nrow = m)
colnames(tpm.matrix) <- datadir

for( i in 1:n ){
    data.i <- read.table(paste0("se_read_mapping/",
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
     file = "se_PCA/output/post_combine_kallisto_results_matrices.RData")


## to rotate matrix: PCA and many other functions require individual in rows and genes in columns
rtpm <- t(apply(tpm.matrix, 2, rev))

rm(list=setdiff(ls(), "rtpm"))


##### ===========================================================================================
##### PCA plotting
##### ===========================================================================================

pdf("se_PCA/output/post_combine_PCA.pdf",height=9,width=9)

pca1=PCA(log2(rtpm+1), graph=F)

datacolor<-as.vector(c("#110B79", "#110B79", "#110B79","#110B79", "#110B79", "#110B79","#110B79", "#110B79", "#110B79","#110B79", "#110B79",
                       "#70AD47", "#70AD47", "#70AD47","#70AD47", "#70AD47", "#70AD47","#70AD47", "#70AD47", "#70AD47",
                       "#70AD47", "#70AD47", "#70AD47","#70AD47", "#70AD47", "#70AD47","#70AD47", "#70AD47"))

plot(pca1$eig[,2], type="p", col="blue", pch=16,
     main="Plot of Principle Components", 
     xlab="Principle Components", 
     ylab="Percentage of Variance")

plot.PCA(pca1, 
         choix="ind",
         ellipse = NULL,
         axes=c(1,2),
         label="none",
         habillage="ind", 
         cex=3,
         title="Plot of Individuals (PC1 vs PC2)",
         palette=palette(datacolor)
)

plot.PCA(pca1, 
         choix="ind",
         axes=c(2,3),
         habillage="ind", 
         label="none", 
         cex=3,
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



## contribution from all variables to PC1 to PC5
dimtank<-as.matrix(pca1$var$contrib)
write.csv(dimtank, 
          file = "se_PCA/output/post_combine_PC1to5_gene_contribution.csv")

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
     file = "se_PCA/output/post_combine_genes_of_PC.RData")

dev.off()

rm(list=ls())
