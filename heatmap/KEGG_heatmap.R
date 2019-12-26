
library(pheatmap)

## working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome") # Set it to your current working directory

## to load existed dataset and skip the dataset generation section
#load("heatmap/output/KEGG_heatmap/UP_matrix.Rdata")
#load("heatmap/output/KEGG_heatmap/DN_matrix.Rdata")

##### ===========================================================================================
##### Dataset generation
##### ===========================================================================================

## to load the results from analyses without P-value cutoff for KEGG (cutoff line P = 1.00), which include  all KEGG even if their P>0.01  
## to keep only 3 columns (KEGG ID, P value, and Adjusted P value) 
D3dUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_D3d_UP1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"), sep ="\t")
D6dUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_D6d_UP1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"),sep ="\t")
D1mUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_D1m_UP1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"),sep ="\t")
NonDUP = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_NonD_UP1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"),sep ="\t")

D3dDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_D3d_DN1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"), sep ="\t")
D6dDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_D6d_DN1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"),sep ="\t")
D1mDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_D1m_DN1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"),sep ="\t")
NonDDN = read.csv("GO_KEGG_enrichment/output/DE_p001/GO_KEGG_p1000/KEGG_NonD_DN1.txt", header = T, colClasses = c("character","numeric","NULL","NULL", "NULL","NULL","NULL","NULL","numeric"),sep ="\t")

## to respectively generate a master UP-regulated KEGG matrix and a master DN-regulated KEGG matrix 
## the UP and DN master matrice will contain the P value of all KEGGs with Padj>0.01 

## to make an index list of 4 conditions (Diapause 3 days, Diapause 6 days, Diapause 1 month, and Non-diapause)
index<-c("D3d", "D6d", "D1m", "NonD")

## doing FDR filtering here.
FDRcutoff<-0.05

## (1) the master UP-regulated matrix
## the obtain a list of UP-regulated KEGGs which are enriched in any of the 4 conditions 
UPKEGG_list<-union(
  union(D3dUP[D3dUP$kpadj<FDRcutoff,1], D6dUP[D6dUP$kpadj<FDRcutoff,1]), 
  union(D1mUP[D1mUP$kpadj<FDRcutoff,1], NonDUP[NonDUP$kpadj<FDRcutoff,1])
)

## to generate a UP-regulated KEGG matrix based on the list
## the input value is -log(P value), with the default value as 0 (P value=1)
UP_matrix<-matrix(data=0, ncol=length(index), nrow=length(UPKEGG_list))
colnames(UP_matrix)<-index
rownames(UP_matrix)<-UPKEGG_list
for (i in 1:length(UPKEGG_list)){
  for (j1 in 1:length(D3dUP[,1]))   {if (UPKEGG_list[i] == D3dUP[j1,1])  {UP_matrix[i,1] <- log(1/D3dUP[j1,2])}}
  for (j2 in 1:length(D6dUP[,1]))   {if (UPKEGG_list[i] == D6dUP[j2,1])  {UP_matrix[i,2] <- log(1/D6dUP[j2,2])}}
  for (j3 in 1:length(D1mUP[,1]))   {if (UPKEGG_list[i] == D1mUP[j3,1])  {UP_matrix[i,3] <- log(1/D1mUP[j3,2])}}
  for (j4 in 1:length(NonDUP[,1]))  {if (UPKEGG_list[i] == NonDUP[j4,1]) {UP_matrix[i,4] <- log(1/NonDUP[j4,2])}}
}

## (2) the master DN-regulated matrix
## the obtain a list of DN-regulated KEGGs which are enriched in any of the 4 conditions 
DNKEGG_list<-union(
  union(D3dDN[D3dDN$kpadj<FDRcutoff,1], D6dDN[D6dDN$kpadj<FDRcutoff,1]), 
  union(D1mDN[D1mDN$kpadj<FDRcutoff,1], NonDDN[NonDDN$kpadj<FDRcutoff,1])
)

## to generate a DN-regulated KEGG matrix based on the list
## the input value is -log(P value), with the default value as 0 (P value=1)
DN_matrix<-matrix(data=0, ncol=length(index), nrow=length(DNKEGG_list))
colnames(DN_matrix)<-index
rownames(DN_matrix)<-DNKEGG_list
for (i in 1:length(DNKEGG_list)){
  for (j5 in 1:length(D3dDN[,1]))   {if (DNKEGG_list[i] == D3dDN[j5,1])  {DN_matrix[i,1] <- log(1/D3dDN[j5,2])}}
  for (j6 in 1:length(D6dDN[,1]))   {if (DNKEGG_list[i] == D6dDN[j6,1])  {DN_matrix[i,2] <- log(1/D6dDN[j6,2])}}
  for (j7 in 1:length(D1mDN[,1]))   {if (DNKEGG_list[i] == D1mDN[j7,1])  {DN_matrix[i,3] <- log(1/D1mDN[j7,2])}}
  for (j8 in 1:length(NonDDN[,1]))  {if (DNKEGG_list[i] == NonDDN[j8,1]) {DN_matrix[i,4] <- log(1/NonDDN[j8,2])}}
}


## (3) to create two lists of KEGG annotation mirroring to the KEGG IDs in both master UP- and DN-regulated KEGG matrice
## to load and match KEGG annotation to both matrice
KEGGlist<-read.csv(file="heatmap/reference/KEGGlist.csv", header=T, colClasses = c("character","character"),sep=',')
UPlist<-as.matrix(KEGGlist[KEGGlist$KEGGID %in% rownames(UP_matrix),])
DNlist<-as.matrix(KEGGlist[KEGGlist$KEGGID %in% rownames(DN_matrix),])

## save to files, raw data to .Rdata, and annotated list to .csv 
save(UP_matrix, UPlist, file = "heatmap/output/KEGG_heatmap/UP_matrix.Rdata")
save(DN_matrix, DNlist, file = "heatmap/output/KEGG_heatmap/DN_matrix.Rdata")
write.csv(UPlist, file = "heatmap/output/KEGG_heatmap/UP_list.csv")
write.csv(DNlist, file = "heatmap/output/KEGG_heatmap/DN_list.csv")


##### ===========================================================================================
##### KEGG annotation appendix 
##### ===========================================================================================

## to create two copies of master UP- and DN-regulated matrice, one using KEGG IDs, and the other KEGG annotations
## to duplicate both master UP- and DN-regulated matrice to two copies 
UPdatamatKEGGID <- as.matrix(UP_matrix) # using KEGG IDs
DNdatamatKEGGID <- as.matrix(DN_matrix) # using KEGG IDs
UPdatamat <- as.matrix(UP_matrix) # to be used for KEGG annotations
DNdatamat <- as.matrix(DN_matrix) # to be used for KEGG annotations

## to reorder the KEGG IDs in UPdatamat & DNdatamat, and then replace KEGG IDs with KEGG annotations
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

## to put a cap on the log (1/P value), setting the dynamic range from 0 to the cap (any P value lower than the cap will be considered the same as the cap) 
## this is to prevent a few extremely low P values to bias the scale and squeeze majority of the KEGGs into a tight window
cutoff=20
UPdatamat[UPdatamat>cutoff]=cutoff
DNdatamat[DNdatamat>cutoff]=cutoff
UPdatamatKEGGID[UPdatamatKEGGID>cutoff]=cutoff
DNdatamatKEGGID[DNdatamatKEGGID>cutoff]=cutoff

color1<-colorRampPalette(c("lightyellow","darkblue")) # for DN-regulated KEGG heatmap
color2<-colorRampPalette(c("lightyellow","red")) # for UP-regulated KEGG heatmap

## to save the results as PDF
pdf(paste0("heatmap/output/KEGG_heatmap/KEGG_overview.pdf"), height=11,width=8.5, paper="letter")
# heatmap of UP-regulated KEGG annotations
pheatmap(UPdatamat, main="UP", Rowv=as.dendrogram(UPclust), clustering_method="ward.D2", scale="none", 
         cluster_cols = F, cellwidth=30, Colv=NA, col=color2(100))
# heatmap of UP-regulated KEGG IDs
pheatmap(UPdatamatKEGGID, main="UP", Rowv=as.dendrogram(UPclust), clustering_method="ward.D2", scale="none", 
         cluster_cols = F, cellwidth=30, Colv=NA, col=color2(100))

# heatmap of DN-regulated KEGG annotations
pheatmap(DNdatamat, main="DN", Rowv=as.dendrogram(DNclust), clustering_method="ward.D2", scale="none",
         cluster_cols = F, cellwidth=30, Colv=NA, col=color1(100))
# heatmap of UP-regulated KEGG IDs
pheatmap(DNdatamatKEGGID, main="DN", Rowv=as.dendrogram(DNclust), clustering_method="ward.D2", scale="none",
         cluster_cols = F, cellwidth=30, Colv=NA, col=color1(100))

dev.off()

