
# Set wd to the current directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

res <- read.csv("se_valcano_plot/BLE_result.csv", header = T, sep = ",")
#head(res)

pdf("se_valcano_plot/valcano.pdf", height=11,width=8.5)


# Make a basic volcano plot
with(res, plot(Log2Fold, -log10(qval), pch=20, main="Volcano plot", xlim=c(-2,2), ylim=c(0,20), col="gray"))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(res, qval<.01 ), points(Log2Fold, -log10(qval), pch=20, col="red"))
#with(subset(res, abs(Log2Fold)>0.58), points(Log2Fold, -log10(qval), pch=20, col="orange"))
#with(subset(res, qval<.01 & abs(Log2Fold)>0.58), points(Log2Fold, -log10(qval), pch=20, col="lightgreen"))

with(subset(res, qval<.01 & Log2Fold>0), points(Log2Fold, -log10(qval), pch=20, col="#FF99CC"))
with(subset(res, qval<.01 & Log2Fold>0.58), points(Log2Fold, -log10(qval), pch=20, col="#FF66B2"))
with(subset(res, qval<.01 & Log2Fold>1), points(Log2Fold, -log10(qval), pch=20, col="#FF0000"))


with(subset(res, qval<.01 & Log2Fold<0), points(Log2Fold, -log10(qval), pch=20, col="#66B2FF"))
with(subset(res, qval<.01 & Log2Fold<(-0.58)), points(Log2Fold, -log10(qval), pch=20, col="#0080FF"))
with(subset(res, qval<.01 & Log2Fold<(-1)), points(Log2Fold, -log10(qval), pch=20, col="#0000FF"))

# Label points with the textxy function from the calibrate plot
#library(calibrate)
#with(subset(res, qval<.001 & abs(Log2Fold)>2), textxy(Log2Fold, -log10(qval), labs=id, cex=.8))


dev.off()

