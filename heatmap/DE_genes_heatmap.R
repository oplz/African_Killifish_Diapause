
##### Codes for the heatmap of all differentially expressed genes 



library(pheatmap)

## to set the "Diapause_transcritome" as the working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

## to load the DE gene list
panelkey<-c("CBX7_DEG_UP_marked_142")
GOlist = read.csv(paste0("gene_lists/DE_heatmap/", panelkey, ".txt"), header = T, col.names=c("id"))

## to load the reference list containing normalized TPMs of all genes
TPMlist = read.csv("heatmap/reference/normalized_TPM_combined.csv", header = T, 
                  colClasses = c("character", 
                                 "numeric","numeric","numeric",
                                 "numeric","numeric","numeric",
                                 "numeric","numeric","numeric",
                                 "numeric","numeric","numeric", 
                                 "numeric","numeric","numeric"), sep =",")
rownames(TPMlist)<-TPMlist$id
TPMlist$id<-NULL

## to make a new list only including the genes on the DE gene list and their TPMs
GOmatrix<-TPMlist[rownames(TPMlist) %in% GOlist$id,]

## to reorder by the gene order from GOlist
#GOmatrix<-GOmatrix[GOlist$id,, drop=F]

## to generate a heatmap
genemat <- as.matrix(GOmatrix)
rownames(genemat)<-rownames(GOmatrix)

GOclust<-hclust(as.dist(1-cor(t(genemat))), method="ward.D2")

colorcode<-colorRampPalette(c("darkblue","lightyellow","red")) 

## to save the results as PDF
pdf(paste0("heatmap/output/", panelkey,".pdf"), height=11,width=8.5, paper="letter")

pheatmap(genemat, main=panelkey, Rowv=as.dendrogram(GOclust), clustering_method="ward.D2", scale="row",
         show_rownames = T, cluster_cols = F, cluster_rows = T, cellwidth=20, fontsize_row=10, Colv=NA, col=colorcode(100))

dev.off()

