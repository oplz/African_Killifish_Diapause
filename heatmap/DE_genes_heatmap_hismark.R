
##### Codes for the heatmap of all differentially expressed genes and histone marking, with the clustering filter



library(pheatmap)

## to set the "Diapause_transcritome" as the working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

## to load the DE gene list
panelkey<-c("--AutophagyNew1_list")
GOlist = read.csv(paste0("gene_lists/DE_heatmap/", panelkey, ".txt"), header = T, col.names=c("id"))
GOlist<-unique(GOlist)

# to load 6 clustering filters
UT = read.csv(paste0("gene_lists/cluster6/UP_throughhout_GeneList_C1.txt"), header = T, col.names=c("id"))
UE = read.csv(paste0("gene_lists/cluster6/UP_early_GeneList_C6.txt"), header = T, col.names=c("id"))
UL = read.csv(paste0("gene_lists/cluster6/UP_late_GeneList_C3.txt"), header = T, col.names=c("id"))
DT = read.csv(paste0("gene_lists/cluster6/DN_throughout_GeneList_C2.txt"), header = T, col.names=c("id"))
DE = read.csv(paste0("gene_lists/cluster6/DN_early_GeneList_C4.txt"), header = T, col.names=c("id"))
DL = read.csv(paste0("gene_lists/cluster6/DN_late_GeneList_C5.txt"), header = T, col.names=c("id"))


## to load the reference list containing normalized TPMs of all genes
DElist = read.csv("heatmap/reference/normalized_TPM_combined.csv", header = T, 
                  colClasses = c("character", 
                                 "numeric","numeric","numeric",
                                 "numeric","numeric","numeric",
                                 "numeric","numeric","numeric",
                                 "numeric","numeric","numeric", 
                                 "numeric","numeric","numeric"), sep =",")
rownames(DElist)<-DElist$id


# to filter the GOlist with preferred clustering filter, keeping only genes matched the clustering pattern
# 6 clustering filter: DT (Down throughout), DE (Down early), DL (Down late), UT (Up throughout), UE (Up early), UL (Up late)
newlist<-as.data.frame(GOlist[GOlist$id %in% UT$id,])


# to obtain RNA-seq data of genes on the newlist
colnames(newlist)<-"assignedorder"
rownames(DElist)<-DElist$id
DElist$id<-NULL
GOmatrix<-DElist[rownames(DElist) %in% newlist$assignedorder,]

#reorder by gene order from newlist
GOmatrix<-GOmatrix[match(newlist$assignedorder,rownames(GOmatrix)),, drop=FALSE]

genemat <- as.matrix(GOmatrix)
rownames(genemat)<-rownames(GOmatrix)


# to load gene list with H3K27me3 histone marks

K27me3DND = read.csv("heatmap/reference/DND_share.txt", header = F, colClasses = c("character"))
K27melost = read.csv("heatmap/reference/lost.txt", header = F, colClasses = c("character"))
K27megained = read.csv("heatmap/reference/gained.txt", header = F, colClasses = c("character"))

# to assign which H3K27me3 list to use
K27list<-K27me3DND
K27list<-K27melost
K27list<-K27megained

# to make a matching list for histone marks, put 0 as devalue value
hismark <- data.frame(data=1, numeric(length = 1 * length(GOmatrix$PreD_Lib1)))
hismark[,1]<-rownames(GOmatrix)
colnames(hismark)<-c("id", "marked_D")
rownames(hismark)<-hismark$id

# to make a matching list of histone mark annotation, if marked, put 10, if not stay as 0.
for (i in 1:length(hismark[,1]))  {if (hismark[i,1] %in% K27list$V1) {hismark[i,2] <- 10}}

GOmatrix$id<-NULL
hismark$id<-NULL

## to generate a heatmap of genes
genemat <- as.matrix(GOmatrix)
rownames(genemat)<-rownames(GOmatrix)
GOclust<-hclust(as.dist(1-cor(t(genemat))), method="ward.D2")
colorcode<-colorRampPalette(c("#0000cc","lightyellow","red")) 

## to do the same to histone marks
genemark <- as.matrix(hismark)
rownames(genemark)<-rownames(hismark)
colorcode2<-colorRampPalette(c("lightyellow","#110b79")) 



## to save the results as PDF
pdf(paste0("heatmap/output/", panelkey,"clustered_lost.pdf"), height=11,width=8.5, paper="letter")

pheatmap(genemat, main=panelkey, Rowv=as.dendrogram(GOclust), clustering_method="ward.D2", scale="row",
         show_rownames = T, cluster_cols = F, cluster_rows = F, cellwidth=20, fontsize_row=10, Colv=NA, col=colorcode(100))

pheatmap(genemark, main=panelkey, Rowv=as.dendrogram(GOclust),
         show_rownames = T, cluster_cols = F, cluster_rows = F, cellwidth=20, fontsize_row=10, Colv=NA, col=colorcode2(100))

dev.off()

# genemark[1,1]<-5 #if makring is universal, position 5 in first row to make it possible for heatmap plotting.  
# This makes the heatmap bar from 5 to 0 or 10, instead of o to 10.
