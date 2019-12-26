library(Rsamtools)
library(parallel)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)
library(biomaRt)
library(ggplot2)
library(pander)
library(plotly)
library(grid)
library(GenomicRanges)
library(reshape)
library(DT)
library(VariantTools)
library(pheatmap)
library(RColorBrewer)
library(ChIPseeker)
library(VennDiagram)
library(genomation)
library(DiffBind)

setwd("/home/wew/Chip_Seq/diapause6/Analysis_H3K4me3")


## Experimental Design
sample_name <- c("D4_H3K4me3_1","D4_H3K4me3_2","ND4_H3K4me3_1","ND4_H3K4me3_2")
input <- c("D4_input","D4_input","ND4_input","ND4_input")
timepoint <- c("diapause","diapause","nondiapause","nondiapause")
antibody <- c("H3K4me3","H3K4me3","H3K4me3","H3K4me3")
replicates <- rep(2, 4)
exp_design <- data.frame(sample_name,input,timepoint,antibody,replicates)
exp_design$order_type <- "ChIP-Seq"
datatable(exp_design)


# Read in the peaksets
bf.dir.1 <- "/home/wew/Chip_Seq/diapause6"
peak.dir.1 <- "/home/wew/Chip_Seq/diapause6/H3K4me3/H3K4me3"
## Make sample sheet
sample.sheet <- data.frame(c("D4_H3K4me3_1","D4_H3K4me3_2","ND4_H3K4me3_1","ND4_H3K4me3_2"), 
                           c("diapause", "diapause","non_diapause", "non_diapause"),
                           c(1,2,1,2),
                           c(paste(bf.dir.1,"H3K4me3_D1_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"H3K4me3_D2_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"H3K4me3_ND1_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"H3K4me3_ND2_sorted_F1540_q15.bam",sep="/")),
                           c("D_input","D_input", "ND_input","ND_input"),
                           c(paste(bf.dir.1,"input_D_K4K27_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"input_D_K4K27_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"input_ND_K4K27_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"input_ND_K4K27_sorted_F1540_q15.bam",sep="/")),
                           c(paste(peak.dir.1,"D1_H3K4me3_macs2_peaks.narrowPeak",sep="/"),
                             paste(peak.dir.1,"D2_H3K4me3_macs2_peaks.narrowPeak",sep="/"),
                             paste(peak.dir.1,"ND1_H3K4me3_macs2_peaks.narrowPeak",sep="/"),
                             paste(peak.dir.1,"ND2_H3K4me3_macs2_peaks.narrowPeak",sep="/")),
                           c(rep("narrow",4)))
colnames(sample.sheet) <- c("SampleID", "Condition", "Replicate", "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller")

## Make DBA object
MyChiP.DBA <- dba(sampleSheet = sample.sheet)
MyChiP.DBA

pdf("H3K4me3_plots.pdf")

plot(MyChiP.DBA)
# Counting reads
MyChiP.DBA.counts <- dba.count(MyChiP.DBA, summits = 200)
plot(MyChiP.DBA.counts)
#Establishing a contrast
MyChiP.DBA.counts.con <- dba.contrast(MyChiP.DBA.counts, categories = DBA_CONDITION, minMembers = 2)
# Performing the differential analysis 
MyChiP.DBA.counts.ana <- dba.analyze(MyChiP.DBA.counts.con)
MyChiP.DBA.counts.ana
dba.plotHeatmap(MyChiP.DBA.counts.ana, contrast = 1) 

dev.off()

# Retrieving the differentially bound sites, I found high p value cut off gives many false positive peaks.
MyChiP.DBA.counts.ana.DB <- dba.report(MyChiP.DBA.counts.ana, th = 0.001, bUsePval = FALSE)
MyChiP.DBA.counts.ana.DB.df <- data.frame(MyChiP.DBA.counts.ana.DB) 

# Retrieving the all bound sites # 44747
MyChiP.DBA.counts.ana.all <- dba.report(MyChiP.DBA.counts.ana, th= 1)
MyChiP.DBA.counts.ana.all.df <- data.frame(MyChiP.DBA.counts.ana.all) 


#* PCA
#* MA plot
#* Box Plot
#* Heatmap

pdf("H3K4me3_PCAMABOX.pdf", 10, 8)
pandoc.header('PCA: all peaks', level=4)
## 4.2 PCA plots
dba.plotPCA(MyChiP.DBA.counts,DBA_CONDITION,label = DBA_CONDITION) # all peaks

pandoc.header('PCA: differential peaks', level=4)
dba.plotPCA(MyChiP.DBA.counts.ana,contrast = 1,label = DBA_CONDITION) # differentially bound sites only 

pandoc.header('MA Plot', level=4)
## 4.3 MA plots
dba.plotMA(MyChiP.DBA.counts.ana)

pandoc.header('Volcano Plot', level=4)
## 4.4 Volcano plots
dba.plotVolcano(MyChiP.DBA.counts.ana)

pandoc.header('BoxPlot Plot', level=4)
## 4.5 Boxplots
sum(MyChiP.DBA.counts.ana.DB$Fold<1)
sum(MyChiP.DBA.counts.ana.DB$Fold>1)
pvals <- dba.plotBox(MyChiP.DBA.counts.ana)

pandoc.header('Heatmaps', level=4)
## 4.6 Heatmaps
corvals <- dba.plotHeatmap(MyChiP.DBA.counts.ana)
corvals <- dba.plotHeatmap(MyChiP.DBA.counts.ana, contrast = 1, correlations = FALSE, scale="row")
dev.off()


##### Annotation

MyChiP.DBA.counts.ana.all 
# annotate peaks
txdb <- loadDb("/home/wew/Nothobranchius_furzeri/dag.analysis/wew2/genome/refseq_r100_plus_wew_08102017_curated_renamed_PPS_CKHU.txdb")
MyChiP.DBA.counts.ana.all.ann <- annotatePeak(MyChiP.DBA.counts.ana.all, TxDb = txdb, tssRegion=c(-3000, 3000))
MyChiP.DBA.counts.ana.all.ann
MyChiP.DBA.counts.ana.all.ann.gr <- as.GRanges(MyChiP.DBA.counts.ana.all.ann)
MyChiP.DBA.counts.ana.all.ann.df <- data.frame(MyChiP.DBA.counts.ana.all.ann.gr) # 69382
#external_order <- BM.info[match(MyChiP.DBA.counts.ana.all.ann.df$geneId, BM.info$ensembl_gene_id),]
#MyChiP.DBA.counts.ana.all.ann.df$external_gene_name <- external_order$external_gene_name
write.table(MyChiP.DBA.counts.ana.all.ann.df,"DiffBind_H3K4me3_all_peaks_annotated.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# filter for peaks up in diapause
MyChiP.DBA.counts.ana.DB.1dpa.up <- MyChiP.DBA.counts.ana.DB[MyChiP.DBA.counts.ana.DB$Fold > 1,]
# annotate peaks
#txdb <- loadDb("/home/wew/Nothobranchius_furzeri/dag.analysis/wew4/genome/danRer10.Ens_91.gff.txdb")
MyChiP.DBA.counts.ana.DB.1dpa.up.ann <- annotatePeak(MyChiP.DBA.counts.ana.DB.1dpa.up, TxDb = txdb, tssRegion=c(-3000, 3000))
MyChiP.DBA.counts.ana.DB.1dpa.up.ann
MyChiP.DBA.counts.ana.DB.1dpa.up.ann.gr <- as.GRanges(MyChiP.DBA.counts.ana.DB.1dpa.up.ann)
MyChiP.DBA.counts.ana.DB.1dpa.up.ann.df <- data.frame(MyChiP.DBA.counts.ana.DB.1dpa.up.ann.gr)
#external_order <- BM.info[match(MyChiP.DBA.counts.ana.DB.1dpa.up.ann.df$geneId, BM.info$ensembl_gene_id),]
#MyChiP.DBA.counts.ana.DB.1dpa.up.ann.df$external_gene_name <- external_order$external_gene_name
write.table(MyChiP.DBA.counts.ana.DB.1dpa.up.ann.df,"DiffBind_H3K4me3_diapause_up.txt", row.names = FALSE, quote = FALSE, sep = "\t")

plotAnnoBar(MyChiP.DBA.counts.ana.DB.1dpa.up.ann)
plotDistToTSS(MyChiP.DBA.counts.ana.DB.1dpa.up.ann)


# filter for peaks down in diapause
MyChiP.DBA.counts.ana.DB.1dpa.down <- MyChiP.DBA.counts.ana.DB[MyChiP.DBA.counts.ana.DB$Fold < -1,]
# annotate peaks
#txdb <- loadDb("/home/wew/Nothobranchius_furzeri/dag.analysis/wew4/genome/danRer10.Ens_91.gff.txdb")
MyChiP.DBA.counts.ana.DB.1dpa.down.ann <- annotatePeak(MyChiP.DBA.counts.ana.DB.1dpa.down, TxDb = txdb, tssRegion=c(-3000, 3000))
MyChiP.DBA.counts.ana.DB.1dpa.down.ann
MyChiP.DBA.counts.ana.DB.1dpa.down.ann.gr <- as.GRanges(MyChiP.DBA.counts.ana.DB.1dpa.down.ann)
MyChiP.DBA.counts.ana.DB.1dpa.down.ann.df <- data.frame(MyChiP.DBA.counts.ana.DB.1dpa.down.ann.gr)
write.table(MyChiP.DBA.counts.ana.DB.1dpa.down.ann.df,"DiffBind_H3K4me3_diapause_down.txt", row.names = FALSE, quote = FALSE, sep = "\t")



##### Heatmap

#diapause enriched peaks

MyChiP.peaks <- resize(MyChiP.DBA.counts.ana.DB.1dpa.up, width = 4000, fix = "center")

# bamfiles 
genomationDataPath.1 <- "/home/wew/Chip_Seq/diapause6"
bam.files.1 <- list.files(genomationDataPath.1, full.names = TRUE, pattern = "bam$")
bam.files <- bam.files.1[grepl("H3K4me3", bam.files.1)]

sml.3 <- ScoreMatrixList(bam.files, MyChiP.peaks, bin.num = 80, type = "bam")
names(sml.3) <- c("D4_H3K4me3_1","D4_H3K4me3_2","ND4_H3K4me3_1","ND4_H3K4me3_1")

pdf("H3K4me3_heatmaptest_diapause_enrichedpeaks.pdf", 10, 8)
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE, 
                col = c("white","purple"))
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE)



plotMeta(sml.3, xcoords = c(-2000, 2000), centralTend = "mean", 
         profile.names = c("D4_H3K4me3_1","D4_H3K4me3_2","ND4_H3K4me3_1","ND4_H3K4me3_1"),
         winsorize=c(0,97))
dev.off()


##### Heatmap

#all peaks

MyChiP.peaks <- resize(MyChiP.DBA.counts.ana.all, width = 4000, fix = "center")
#MyChiP.peaks <- MyChiP.DBA.counts.ana.DB
  
# bamfiles 
genomationDataPath.1 <- "/home/wew/Chip_Seq/diapause6"
bam.files.1 <- list.files(genomationDataPath.1, full.names = TRUE, pattern = "bam$")
bam.files <- bam.files.1[grepl("H3K4me3", bam.files.1)]

sml.3 <- ScoreMatrixList(bam.files, MyChiP.peaks, bin.num = 80, type = "bam")
names(sml.3) <- c("D4_H3K4me3_1","D4_H3K4me3_2","ND4_H3K4me3_1","ND4_H3K4me3_1")

pdf("H3K4me3_heatmaptest_diapause_alldf_peaks.pdf", 10, 8)
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE, 
                col = c("white", "pink"))
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE)



plotMeta(sml.3, xcoords = c(-2000, 2000), centralTend = "mean", 
         profile.names = c("D4_H3K4me3_1","D4_H3K4me3_2","ND4_H3K4me3_1","ND4_H3K4me3_1"),
         winsorize=c(0,97))

dev.off()
