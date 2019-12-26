
library(pheatmap)

## working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

## to load existed dataset and skip the dataset generation section
#load("heatmap/output/GO_term_heatmap/UP_matrix.Rdata")
#load("heatmap/output/GO_term_heatmap/DN_matrix.Rdata")

##### ===========================================================================================
##### Dataset generation
##### ===========================================================================================

## to load the results from analyses without P-value cutoff for GO terms (cutoff line P = 1.00), which include  all GO terms even if their P>0.01  
## to keep only 3 columns (GO term ID, P value, and Adjusted P value) 
D3dUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_D3d_UP1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"), sep ="\t")
D6dUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_D6d_UP1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
D1mUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_D1m_UP1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
NonDUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_NonD_UP1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")

D3dDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_D3d_DN1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"), sep ="\t")
D6dDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_D6d_DN1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
D1mDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_D1m_DN1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")
NonDDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/GO_NonD_DN1.txt", header = T, colClasses = c("character","numeric","NULL", "NULL","NULL","NULL","NULL","NULL","numeric"),sep ="\t")


## to respectively generate a master UP-regulated GO term matrix and a master DN-regulated GO term matrix 
## the UP and DN master matrice will contain the P value of all GO terms with Padj>0.01 

## to make an index list of 4 conditions (Diapause 3 days, Diapause 6 days, Diapause 1 month, and Non-diapause)
index<-c("D3d", "D6d", "D1m", "NonD")

## doing FDR filtering here.
FDRcutoff<-0.01


## (1) the master UP-regulated matrix
## the obtain a list of UP-regulated GO terms which are enriched in any of the 4 conditions 
UPGOterm_list<-union(
  union(D3dUP[D3dUP$padj<FDRcutoff,1], D6dUP[D6dUP$padj<FDRcutoff,1]), 
  union(D1mUP[D1mUP$padj<FDRcutoff,1], NonDUP[NonDUP$padj<FDRcutoff,1])
)

## to generate a UP-regulated GO term matrix based on the list
## the input value is -log(P value), with the default value as 0 (P value=1)
UP_matrix<-matrix(data=0, ncol=length(index), nrow=length(UPGOterm_list))
colnames(UP_matrix)<-index
rownames(UP_matrix)<-UPGOterm_list
for (i in 1:length(UPGOterm_list)){
  for (j1 in 1:length(D3dUP[,1]))   {if (UPGOterm_list[i] == D3dUP[j1,1])  {UP_matrix[i,1] <- log(1/D3dUP[j1,2])}}
  for (j2 in 1:length(D6dUP[,1]))   {if (UPGOterm_list[i] == D6dUP[j2,1])  {UP_matrix[i,2] <- log(1/D6dUP[j2,2])}}
  for (j3 in 1:length(D1mUP[,1]))   {if (UPGOterm_list[i] == D1mUP[j3,1])  {UP_matrix[i,3] <- log(1/D1mUP[j3,2])}}
  for (j4 in 1:length(NonDUP[,1]))  {if (UPGOterm_list[i] == NonDUP[j4,1]) {UP_matrix[i,4] <- log(1/NonDUP[j4,2])}}
}
  
## (2) the master DN-regulated matrix
## the obtain a list of DN-regulated GO terms which are enriched in any of the 4 conditions 
DNGOterm_list<-union(
  union(D3dDN[D3dDN$padj<FDRcutoff,1], D6dDN[D6dDN$padj<FDRcutoff,1]), 
  union(D1mDN[D1mDN$padj<FDRcutoff,1], NonDDN[NonDDN$padj<FDRcutoff,1])
)

## to generate a DN-regulated GO term matrix based on the list
## the input value is -log(P value), with the default value as 0 (P value=1)
DN_matrix<-matrix(data=0, ncol=length(index), nrow=length(DNGOterm_list))
colnames(DN_matrix)<-index
rownames(DN_matrix)<-DNGOterm_list
for (i in 1:length(DNGOterm_list)){
  for (j5 in 1:length(D3dDN[,1]))   {if (DNGOterm_list[i] == D3dDN[j5,1])  {DN_matrix[i,1] <- log(1/D3dDN[j5,2])}}
  for (j6 in 1:length(D6dDN[,1]))   {if (DNGOterm_list[i] == D6dDN[j6,1])  {DN_matrix[i,2] <- log(1/D6dDN[j6,2])}}
  for (j7 in 1:length(D1mDN[,1]))   {if (DNGOterm_list[i] == D1mDN[j7,1])  {DN_matrix[i,3] <- log(1/D1mDN[j7,2])}}
  for (j8 in 1:length(NonDDN[,1]))  {if (DNGOterm_list[i] == NonDDN[j8,1]) {DN_matrix[i,4] <- log(1/NonDDN[j8,2])}}
}


## (3) to create two lists of GO term annotation mirroring to the GO term IDs in both master UP- and DN-regulated GO term matrice
## to load and match GO term annotation to both matrice
GOlist<-read.csv(file="heatmap/reference/GOlist.csv", header=T, sep=',')
UPlist<-as.matrix(GOlist[GOlist$GOBPID %in% rownames(UP_matrix),])
DNlist<-as.matrix(GOlist[GOlist$GOBPID %in% rownames(DN_matrix),])

## to remove redundants if any
UPlist<-UPlist[!duplicated(UPlist[,1]),]
DNlist<-DNlist[!duplicated(DNlist[,1]),]


## save to files, raw data to .Rdata, and annotated list to .csv 
save(UP_matrix, UPlist, file = "heatmap/output/GO_term_heatmap/UP_matrix.Rdata")
save(DN_matrix, DNlist, file = "heatmap/output/GO_term_heatmap/DN_matrix.Rdata")
write.csv(UPlist, file = "heatmap/output/GO_term_heatmap/UP_list.csv")
write.csv(DNlist, file = "heatmap/output/GO_term_heatmap/DN_list.csv")


##### ===========================================================================================
##### GO term annotation appendix 
##### ===========================================================================================

## to create two copies of master UP- and DN-regulated matrice, one using GO term IDs, and the other GO term annotations
## to duplicate both master UP- and DN-regulated matrice to two copies 
UPdatamatGOID <- as.matrix(UP_matrix) # using GO term IDs
DNdatamatGOID <- as.matrix(DN_matrix) # using GO term IDs
UPdatamat <- as.matrix(UP_matrix) # to be used for GO term annotation
DNdatamat <- as.matrix(DN_matrix) # to be used for GO term annotation

## to reorder the GO term IDs in UPdatamat & DNdatamat, and then replace GO term IDs with GO term annotations
UPdatamat<-UPdatamat[UPlist[,1],, drop=F]
DNdatamat<-DNdatamat[DNlist[,1],, drop=F]
rownames(UPdatamat)<-UPlist[,2]
rownames(DNdatamat)<-DNlist[,2]



##### ===========================================================================================
##### Heatmap generation
##### ===========================================================================================

## clustering
UPclust<-hclust(as.dist(1-cor(t(UPdatamat))), method="ward.D2")
DNclust<-hclust(as.dist(1-cor(t(DNdatamat))), method="ward.D2")

#write.csv(UPdatamat, file = "heatmap/output/GO_term_heatmap/UPdatamatrix.csv")
#write.csv(DNdatamat, file = "heatmap/output/GO_term_heatmap/DNdatamatrix.csv")

## to put a cap on the log (1/P value), setting the dynamic range from 0 to the cap (any P value lower than the cap will be considered the same as the cap) 
## this is to prevent a few extremely low P values to bias the scale and squeeze majority of the GO terms into a tight window
cutoff=25 
UPdatamat[UPdatamat>cutoff]=cutoff
DNdatamat[DNdatamat>cutoff]=cutoff
UPdatamatGOID[UPdatamatGOID>cutoff]=cutoff
DNdatamatGOID[DNdatamatGOID>cutoff]=cutoff

color1<-colorRampPalette(c("lightyellow","darkblue")) # for DN-regulated GO term heatmap
color2<-colorRampPalette(c("lightyellow","red")) # for UP-regulated GO term heatmap

## to save the results as PDF
pdf(paste0("heatmap/output/GO_term_heatmap/GO_overview.pdf"), height=11,width=8.5, paper="letter")

# heatmap of UP-regulated GO term annotations
pheatmap(UPdatamat, main="UP", Rowv=as.dendrogram(UPclust), clustering_method="ward.D2", scale="none", 
         cluster_cols = F, cellwidth=30, fontsize_row=2, Colv=NA, col=color2(100))
# heatmap of UP-regulated GO term IDs
pheatmap(UPdatamatGOID, main="UP", Rowv=as.dendrogram(UPclust), clustering_method="ward.D2", scale="none", 
         cluster_cols = F, cellwidth=30, fontsize_row=2, Colv=NA, col=color2(100))

# heatmap of UP-regulated GO term annotations
pheatmap(DNdatamat, main="DN", Rowv=as.dendrogram(DNclust), clustering_method="ward.D2", scale="none",
         cluster_cols = F, cellwidth=30, fontsize_row=1, Colv=NA, col=color1(100))
# heatmap of DN-regulated GO term IDs
pheatmap(DNdatamatGOID, main="DN", Rowv=as.dendrogram(DNclust), clustering_method="ward.D2", scale="none",
         cluster_cols = F, cellwidth=30, fontsize_row=1, Colv=NA, col=color1(100))

dev.off()

