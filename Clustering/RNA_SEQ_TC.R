###-------------------------------------------------------------------------------------------------------------###
####Used to generate clusters of RNAseq data###
###Make sure to 1)set working directory, 2)set input file, and 3)set output file###

#Setup Enviroment
library(TCseq)
###(1)Working Directory###
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome/Clustering")   # Set it to your current working directory
set.seed(1234)

#setup input Data
Ddesign <- read.table("Input/diapause_design_real.txt", header=TRUE)

###(2)Select which data set to run by unhashing it###
#All data
Dfeatures <- read.table("Input/diapause_features_real.txt", header=TRUE)
Dcounts <- read.table("Input/diapause_counts_real.txt", header=TRUE, row.names = 1)

Dcounts <- data.matrix(Dcounts)

#Generate Timecourse Object
Dtca <- TCA(design=Ddesign, genomicFeature=Dfeatures, counts=Dcounts)
tcaDA <- DBanalysis(Dtca, filter.type=NULL)
tcaTT <- timecourseTable(tcaDA, value="FC") #can subset here
tcaTC <- timeclust(tcaTT, algo="cm", k=6, standardize=TRUE) #change the value of K to change cluster number

#Plot
tcaPT <- timeclustplot(tcaTC, value="z-score(Log2FC)", col=3)

#Extract data from timecourse object
ClusterT <- tcaTC@clusterRes@cluster
ZscoreT <- tcaTC@clusterRes@data

###(3)Select appropriate output file by unhashing it###
#All Data
write.table(ClusterT, file="Output/Gene_Cluster6_real.txt", sep="\t")
write.table(ZscoreT, file="Output/Gene_Zscore6_real.txt", sep="\t")

