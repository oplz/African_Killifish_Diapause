# to plot Manhanttan plot of DEGs
library(ggplot2)

setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

df = as.data.frame(read.csv("DE_Manhattan_plot/data.csv", header = T))

# Adjust colors
DE_colors <- c("#66b2ff", "#0080ff", "#004c99", "#ff9933")
names(DE_colors) <- levels(df$Condition)

# Make custom color column to facilitate grey coloring by threshold.
col <- DE_colors[df$Condition]
col[df$FDR > 0.01] <- "#D3D3D3"
df$col <- as.factor(col)

q <- ggplot(df, aes(x = Condition, y = DE, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Differential Expression Folds") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "DE folds (Log2)") +
  theme_minimal() +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")

pdf(paste0("DE_Manhattan_plot/output/DEGs_Manhattan_plot.pdf"), height=11,width=8.5, paper="letter")

q

dev.off()
