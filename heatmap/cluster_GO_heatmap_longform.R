
library(pheatmap)

## working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

##### ===========================================================================================
##### Dataset generation
##### ===========================================================================================

## to load the results from analyses without P-value cutoff for GO terms (cutoff line P = 1.00), which include  all GO terms even if their P>0.01  
## to keep only 3 columns (GO term ID, P value, and Adjusted P value) 
UP_T = read.csv("GO_KEGG_enrichment/output/cluster6_p100/GO_GeneList_C1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"), sep ="\t")
UP_E = read.csv("GO_KEGG_enrichment/output/cluster6_p100/GO_GeneList_C6.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
UP_L = read.csv("GO_KEGG_enrichment/output/cluster6_p100/GO_GeneList_C3.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
DN_T = read.csv("GO_KEGG_enrichment/output/cluster6_p100/GO_GeneList_C2.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"), sep ="\t")
DN_E = read.csv("GO_KEGG_enrichment/output/cluster6_p100/GO_GeneList_C4.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
DN_L = read.csv("GO_KEGG_enrichment/output/cluster6_p100/GO_GeneList_C5.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")


## to respectively generate a master UP-regulated GO term matrix and a master DN-regulated GO term matrix 
## the UP and DN master matrice will contain the P value of all GO terms with Padj>0.01 

## to make an index list of 4 conditions (Diapause 3 days, Diapause 6 days, Diapause 1 month, and Non-diapause)
index<-c("UP_Throughout", "UP_early", "UP_late", "DN_Throughout", "DN_early", "DN_late")

## doing FDR filtering here.
FDRcutoff<-0.001


## (1) the master matrix
## the obtain a list of terms which are enriched 
UPGOterm_list<-union(
  UP_T[UP_T$padj<FDRcutoff,1],
  union(UP_E[UP_E$padj<FDRcutoff,1], UP_L[UP_L$padj<FDRcutoff,1])
  )

DNGOterm_list<-union(
  DN_T[DN_T$padj<FDRcutoff,1],
  union(DN_E[UP_E$padj<FDRcutoff,1], DN_L[DN_L$padj<FDRcutoff,1])
)

GOterm_list<-union(UPGOterm_list, DNGOterm_list)


## to generate a UP-regulated GO term matrix based on the list
## the input value is -log(P value), with the default value as 0 (P value=1)
GO_matrix<-matrix(data=0, ncol=length(index), nrow=length(GOterm_list))
colnames(GO_matrix)<-index
rownames(GO_matrix)<-GOterm_list
for (i in 1:length(GOterm_list)){
  for (j1 in 1:length(UP_T[,1]))   {if (GOterm_list[i] == UP_T[j1,1])  {GO_matrix[i,1] <- log(1/UP_T[j1,2])}}
  for (j2 in 1:length(UP_E[,1]))   {if (GOterm_list[i] == UP_E[j2,1])  {GO_matrix[i,2] <- log(1/UP_E[j2,2])}}
  for (j3 in 1:length(UP_L[,1]))   {if (GOterm_list[i] == UP_L[j3,1])  {GO_matrix[i,3] <- log(1/UP_L[j3,2])}}
  for (j4 in 1:length(DN_T[,1]))   {if (GOterm_list[i] == DN_T[j4,1])  {GO_matrix[i,4] <- log(1/DN_T[j4,2])}}
  for (j5 in 1:length(DN_E[,1]))   {if (GOterm_list[i] == DN_E[j5,1])  {GO_matrix[i,5] <- log(1/DN_E[j5,2])}}
  for (j6 in 1:length(DN_L[,1]))   {if (GOterm_list[i] == DN_L[j6,1])  {GO_matrix[i,6] <- log(1/DN_L[j6,2])}}
}
  

## (3) to create two lists of GO term annotation mirroring to the GO term IDs in both master UP- and DN-regulated GO term matrice
## to load and match GO term annotation to both matrice
GOlist<-read.csv(file="heatmap/reference/GOlist.csv", header=T, sep=',')
cluster_GOlist<-as.matrix(GOlist[GOlist$GOBPID %in% rownames(GO_matrix),])

## to remove redundants if any
cluster_GOlist<-cluster_GOlist[!duplicated(cluster_GOlist[,1]),]



## save to files, raw data to .Rdata, and annotated list to .csv 
save(GO_matrix, cluster_GOlist, file = "heatmap/output/cluster_GO_term_heatmap/cluster_GO_matrix.Rdata")
write.csv(cluster_GOlist, file = "heatmap/output/cluster_GO_term_heatmap/cluster_GOlist.csv")



##### ===========================================================================================
##### GO term annotation appendix 
##### ===========================================================================================

## to create two copies of master matrice, one using GO term IDs, and the other GO term annotations
GOdatamatGOID <- as.matrix(GO_matrix) # using GO term IDs
GOdatamat <- as.matrix(GO_matrix) # to be used for GO term annotation

## to reorder the GO term IDs, and then replace GO term IDs with GO term annotations
GOdatamat<-GOdatamat[cluster_GOlist[,1],, drop=F]
rownames(GOdatamat)<-cluster_GOlist[,2]




##### ===========================================================================================
##### Heatmap generation
##### ===========================================================================================

color1<-colorRampPalette(c("lightyellow","darkblue")) # for DN-regulated GO term heatmap
color2<-colorRampPalette(c("lightyellow","red")) # for UP-regulated GO term heatmap

## clustering
GOclust<-hclust(as.dist(1-cor(t(GOdatamat))), method="ward.D2")

#write.csv(GOdatamat, file = "heatmap/output/cluster_GO_term_heatmap/cluster_GOdatamatrix.csv")

## to put a cap on the log (1/P value), setting the dynamic range from 0 to the cap (any P value lower than the cap will be considered the same as the cap) 
## this is to prevent a few extremely low P values to bias the scale and squeeze majority of the GO terms into a tight window
cutoff=25
GOdatamat[GOdatamat>cutoff]=cutoff
GOdatamatGOID[GOdatamatGOID>cutoff]=cutoff


## to save the results as PDF
pdf(paste0("heatmap/output/cluster_GO_term_heatmap/cluster_GO_overview.pdf"), height=40,width=8.5)

# heatmap of GO term annotations
pheatmap(GOdatamat, main="Cluster 6 GO terms", Rowv=as.dendrogram(GOclust), clustering_method="ward.D2", scale="none", 
         cluster_cols = F, cellwidth=30, fontsize_row=2, Colv=NA, col=color2(100))
# heatmap of GO term IDs
pheatmap(GOdatamatGOID, main="Cluster 6 GO term IDs", Rowv=as.dendrogram(GOclust), clustering_method="ward.D2", scale="none", 
         cluster_cols = F, cellwidth=30, fontsize_row=2, Colv=NA, col=color2(100))

dev.off()

