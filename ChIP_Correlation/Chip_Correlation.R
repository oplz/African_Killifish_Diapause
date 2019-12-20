
# To setup enviroment and load files
library(Hmisc)
set.seed(1234)
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome") # Set it to your current working directory

data1 <- read.table("ChIP_Correlation/reference/k4_vs_rna.txt", header=T, stringsAsFactors = FALSE)
data2 <- read.table("ChIP_Correlation/reference/k27_vs_rna.txt", header=T, stringsAsFactors = FALSE)
data4 <- read.table("ChIP_Correlation/reference/sc_vs_bulk_double_DEG.txt", header=T, stringsAsFactors = FALSE)

# To choose which data set to use 
current <- data2


## OPTIONAL, plotting only a subset of genes 
## To clean up 'current' and keep only interested genes from another list (data3)
  # The panel name for input
  panelkey<-c("LK27_UPE134")

  data3 <- read.table(paste0("gene_lists/ChIPseq/", panelkey, ".txt"), header=T, stringsAsFactors = FALSE)

  # To make an index list of genes present in both current and data3
  index<-as.matrix(intersect(current$Gene, data3$id))

  # To make data frame of interested gene from index
  UP_matrix<-data.frame(data=0, ncol=length(current[1,]), nrow=length(index))
  colnames(UP_matrix)<-colnames(current)

  for (i in 1:length(index)){
    for (j1 in 1:length(current[,1]))   
    {if (index[i] == current[j1,1])  
      {UP_matrix[i,1] <- current[j1,1]
      UP_matrix[i,2] <- current[j1,2]
      UP_matrix[i,3] <- current[j1,3]
      }}
  }
  
  # update current
  current<-UP_matrix




#Create Variables
current$Color <- NA
outx <- 5
outy <- 6
Color <- 4
names <- as.character(current$Gene)
x <- current$Chip_FC
y <- current$RNA_FC

current[,Color] <- "grey"

#label outliers
xers <- as.character(outliers::scores(x, type="chisq", prob=0.95))
yers <- as.character(outliers::scores(y, type="chisq", prob=0.95))
current <- cbind(current,xers,yers)

cor.test(x,y)
cor.test(x,y, method = 'spearman')
summary(lm(x ~ y))

par(mar =c(6,6,7,2))
plot(y,x, main = "Correlation clustered genes", xlab="bulkRNA Log2FC", ylab="ChIP-seq H3K27me3 Log2FC", col=current$Color, pch=19, xlim = c(-8,8), ylim=c(-8,8), cex.main=2, cex.lab=2, cex.axis=1.75)
abline(h = 0)
abline(v = 0)
abline(lm(x ~ y), col="Black", lwd=5)


