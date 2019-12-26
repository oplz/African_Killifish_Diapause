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
sample_name <- c("ND4_H3K4me3_1","ND4_H3K4me3_2","D4_H3K4me3_1","D4_H3K4me3_2")
input <- c("ND4_input","ND4_input", "D4_input","D4_input")
timepoint <- c("nondiapause","nondiapause","diapause","diapause")
antibody <- c("H3K4me3","H3K4me3","H3K4me3","H3K4me3")
replicates <- rep(2, 4)
exp_design <- data.frame(sample_name,input,timepoint,antibody,replicates)
exp_design$order_type <- "ChIP-Seq"
datatable(exp_design)


# Read in the peaksets
bf.dir.1 <- "/home/wew/Chip_Seq/diapause6"
peak.dir.1 <- "/home/wew/Chip_Seq/diapause6/H3K4me3/H3K4me3"
## Make sample sheet
sample.sheet <- data.frame(c("ND4_H3K4me3_1","ND4_H3K4me3_2","D4_H3K4me3_1","D4_H3K4me3_2"), 
                           c("non_diapause", "non_diapause","diapause", "diapause"),
                           c(1,2,1,2),
                           c(paste(bf.dir.1,"H3K4me3_ND1_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"H3K4me3_ND2_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"H3K4me3_D1_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"H3K4me3_D2_sorted_F1540_q15.bam",sep="/")),
                           c("ND4_input","ND4_input", "D4_input","D4_input"),
                           c(paste(bf.dir.1,"input_ND_K4K27_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"input_ND_K4K27_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"input_D_K4K27_sorted_F1540_q15.bam",sep="/"),
                             paste(bf.dir.1,"input_D_K4K27_sorted_F1540_q15.bam",sep="/")),
                           c(paste(peak.dir.1,"ND1_H3K4me3_macs2_peaks.narrowPeak",sep="/"),
                             paste(peak.dir.1,"ND2_H3K4me3_macs2_peaks.narrowPeak",sep="/"),
                             paste(peak.dir.1,"D1_H3K4me3_macs2_peaks.narrowPeak",sep="/"),
                             paste(peak.dir.1,"D2_H3K4me3_macs2_peaks.narrowPeak",sep="/")),
                           c(rep("narrow",4)))
colnames(sample.sheet) <- c("SampleID", "Condition", "Replicate", "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller")

## Make DBA object
MyChIP.DBA <- dba(sampleSheet = sample.sheet)
MyChIP.DBA

pdf("H3K4me3_plotsa.pdf")

plot(MyChIP.DBA)

# Counting reads
MyChIP.DBA.counts <- dba.count(MyChIP.DBA, summits = 400)

plot(MyChIP.DBA.counts)

# Establishing a contrast
MyChIP.DBA.counts.con <- dba.contrast(MyChIP.DBA.counts, categories = DBA_CONDITION, minMembers = 2)

# Performing the differential analysis 
MyChIP.DBA.counts.ana <- dba.analyze(MyChIP.DBA.counts.con)
MyChIP.DBA.counts.ana
dba.plotHeatmap(MyChIP.DBA.counts.ana, contrast = 1) 

dev.off()

# Retrieving the differentially bound sites
MyChIP.DBA.counts.ana.DB <- dba.report(MyChIP.DBA.counts.ana, th = 0.001, bUsePval = FALSE)
MyChIP.DBA.counts.ana.DB.df <- data.frame(MyChIP.DBA.counts.ana.DB) 

# Retrieving the all bound sites # 44747
MyChIP.DBA.counts.ana.all <- dba.report(MyChIP.DBA.counts.ana, th= 1)
MyChIP.DBA.counts.ana.all.df <- data.frame(MyChIP.DBA.counts.ana.all) 


#* PCA: Larger difference between 1dpa and 0dpa when looking at just differential peaks
#* MA plot: pink- Log FC of at least 2 and FDR < 0.01
#* Box Plot: shows that there is a larger difference between 0dpa and 1dpa when peaks are up in 1 dpa
#* Heatmap: small cluster of peaks that are up in 1dpa

pdf("H3K4me3_allplots1.pdf", 10, 8)
pandoc.header('PCA: all peaks', level=4)
## 4.2 PCA plots
dba.plotPCA(MyChIP.DBA.counts,DBA_CONDITION,label = DBA_CONDITION) # all peaks

pandoc.header('PCA: differential peaks', level=4)
dba.plotPCA(MyChIP.DBA.counts.ana,contrast = 1,label = DBA_CONDITION) # differentially bound sites only 

pandoc.header('MA Plot', level=4)
## 4.3 MA plots
dba.plotMA(MyChIP.DBA.counts.ana)

pandoc.header('Volcano Plot', level=4)
## 4.4 Volcano plots
dba.plotVolcano(MyChIP.DBA.counts.ana)

pandoc.header('BoxPlot Plot', level=4)
## 4.5 Boxplots
sum(MyChIP.DBA.counts.ana.DB$Fold<1)
sum(MyChIP.DBA.counts.ana.DB$Fold>1)
pvals <- dba.plotBox(MyChIP.DBA.counts.ana)

pandoc.header('Heatmaps', level=4)
## 4.6 Heatmaps
corvals <- dba.plotHeatmap(MyChIP.DBA.counts.ana)
corvals <- dba.plotHeatmap(MyChIP.DBA.counts.ana, contrast = 1, correlations = FALSE, scale="row")
dev.off()


#### All peaks Up in 1dpa 

#* Peaks up in 1dpa:
#    + all: 12,575
#    + FDR < 0.01: 5831
#* annotation saved: /home/wew/Nothobranchius_furzeri/dag.analysis/wew4/results/annotated_peaks/annotated_DiffBind/

##### Annotation

#MyChIP.DBA.counts.ana.all.notsig <- MyChIP.DBA.counts.ana.all[MyChIP.DBA.counts.ana.all$FDR > 0.01,] # 38867
# annotate peaks
#txdb <- loadDb("/home/wew/Nothobranchius_furzeri/dag.analysis/wew2/genome/refseq_r100_plus_wew_08102017_curated_renamed_PPS_CKHU.txdb")

#MyChIP.DBA.counts.ana.all.notsig.ann <- annotatePeak(MyChIP.DBA.counts.ana.all.notsig, TxDb = txdb, tssRegion=c(-3000, 3000))
#MyChIP.DBA.counts.ana.all.notsig.ann
#MyChIP.DBA.counts.ana.all.notsig.ann.gr <- as.GRanges(MyChIP.DBA.counts.ana.all.notsig.ann)
#MyChIP.DBA.counts.ana.all.notsig.ann.df <- data.frame(MyChIP.DBA.counts.ana.all.notsig.ann.gr) # 38671
#external_order <- BM.info[match(MyChIP.DBA.counts.ana.all.notsig.ann.df$geneId, BM.info$ensembl_gene_id),]
#MyChIP.DBA.counts.ana.all.notsig.ann.df$external_gene_name <- external_order$external_gene_name
#write.table(MyChIP.DBA.counts.ana.all.notsig.ann.df,"DiffBind_H3K4me3_greater_fdr01_notsig.txt", row.names = FALSE, quote = FALSE, sep = "\t")


MyChIP.DBA.counts.ana.all # 69699
# annotate peaks
txdb <- loadDb("/home/wew/Nothobranchius_furzeri/dag.analysis/wew2/genome/refseq_r100_plus_wew_08102017_curated_renamed_PPS_CKHU.txdb")
MyChIP.DBA.counts.ana.all.ann <- annotatePeak(MyChIP.DBA.counts.ana.all, TxDb = txdb, tssRegion=c(-3000, 3000))
MyChIP.DBA.counts.ana.all.ann
MyChIP.DBA.counts.ana.all.ann.gr <- as.GRanges(MyChIP.DBA.counts.ana.all.ann)
MyChIP.DBA.counts.ana.all.ann.df <- data.frame(MyChIP.DBA.counts.ana.all.ann.gr) # 69382
#external_order <- BM.info[match(MyChIP.DBA.counts.ana.all.ann.df$geneId, BM.info$ensembl_gene_id),]
#MyChIP.DBA.counts.ana.all.ann.df$external_gene_name <- external_order$external_gene_name
write.table(MyChIP.DBA.counts.ana.all.ann.df,"DiffBind_H3K4me3_all_peaks_nondiapause.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# filter for peaks up in 1dpa
MyChIP.DBA.counts.ana.DB.1dpa.up <- MyChIP.DBA.counts.ana.DB[MyChIP.DBA.counts.ana.DB$Fold > 1,]
# annotate peaks
#txdb <- loadDb("/home/wew/Nothobranchius_furzeri/dag.analysis/wew4/genome/danRer10.Ens_91.gff.txdb")
MyChIP.DBA.counts.ana.DB.1dpa.up.ann <- annotatePeak(MyChIP.DBA.counts.ana.DB.1dpa.up, TxDb = txdb, tssRegion=c(-3000, 3000))
MyChIP.DBA.counts.ana.DB.1dpa.up.ann
MyChIP.DBA.counts.ana.DB.1dpa.up.ann.gr <- as.GRanges(MyChIP.DBA.counts.ana.DB.1dpa.up.ann)
MyChIP.DBA.counts.ana.DB.1dpa.up.ann.df <- data.frame(MyChIP.DBA.counts.ana.DB.1dpa.up.ann.gr)
#external_order <- BM.info[match(MyChIP.DBA.counts.ana.DB.1dpa.up.ann.df$geneId, BM.info$ensembl_gene_id),]
#MyChIP.DBA.counts.ana.DB.1dpa.up.ann.df$external_gene_name <- external_order$external_gene_name
write.table(MyChIP.DBA.counts.ana.DB.1dpa.up.ann.df,"DiffBind_H3K4me3_fdr01_up_nondiapause.txt", row.names = FALSE, quote = FALSE, sep = "\t")

plotAnnoBar(MyChIP.DBA.counts.ana.DB.1dpa.up.ann)

plotDistToTSS(MyChIP.DBA.counts.ana.DB.1dpa.up.ann)


##### Heatmap

#* Heatmap looking at coverage +/- 2 kb around centered peak
#* Bin size = 80

MyChIP.peaks <- resize(MyChIP.DBA.counts.ana.DB.1dpa.up, width = 4000, fix = "center")
#MyChIP.peaks <- MyChIP.DBA.counts.ana.DB.1dpa.up.FC

# bamfiles 
genomationDataPath.1 <- "/home/wew/Chip_Seq/diapause6"
bam.files.1 <- list.files(genomationDataPath.1, full.names = TRUE, pattern = "bam$")
bam.files <- bam.files.1[grepl("H3K4me3", bam.files.1)]

sml.3 <- ScoreMatrixList(bam.files, MyChIP.peaks, bin.num = 80, type = "bam")
names(sml.3) <- c("D4_H3K4me3_1","D4_H3K4me3_2", "ND4_H3K4me3_1","ND4_H3K4me3_1")

pdf("H3K4me3_heatmap_nondiapause.pdf", 10, 8)
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE, 
                col = c("white", "pink"))
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE)



plotMeta(sml.3, xcoords = c(-2000, 2000), centralTend = "mean", 
         profile.names = c("D4_H3K4me3_1","D4_H3K4me3_2", "ND4_H3K4me3_1","ND4_H3K4me3_1"),
         winsorize=c(0,97))
dev.off()
######################################

MyChIP.peaks <- resize(MyChIP.DBA.counts.ana.all, width = 4000, fix = "center")
#MyChIP.peaks <- MyChIP.DBA.counts.ana.DB

# bamfiles 
genomationDataPath.1 <- "/home/wew/Chip_Seq/diapause6"
bam.files.1 <- list.files(genomationDataPath.1, full.names = TRUE, pattern = "bam$")
bam.files <- bam.files.1[grepl("H3K4me3", bam.files.1)]

sml.3 <- ScoreMatrixList(bam.files, MyChIP.peaks, bin.num = 80, type = "bam")
names(sml.3) <- c("D4_H3K4me3_1","D4_H3K4me3_2", "ND4_H3K4me3_1","ND4_H3K4me3_1")

pdf("H3K4me3_heatmaptest_global_peaks.pdf", 10, 8)
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE, 
                col = c("white", "pink"))
multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE)



plotMeta(sml.3, xcoords = c(-2000, 2000), centralTend = "mean", 
         profile.names = c("D4_H3K4me3_1","D4_H3K4me3_2", "ND4_H3K4me3_1","ND4_H3K4me3_1"),
         winsorize=c(0,97))

dev.off()
