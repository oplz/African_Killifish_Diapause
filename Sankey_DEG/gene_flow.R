## Sankey diagram to track dynamic of DEG compisition throughout conditions

library(ggalluvial)
library(ggplot2)

setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")  # Set it to your current working directory

redtable<-data.frame(read.csv("Sankey_DEG/DEG_reference.csv"))

#vignette(topic = "ggalluvial", package = "ggalluvial")

# releveling the order
redtable$b_D3d <- factor(redtable$b_D3d, levels = c("UP", "NonDE", "DN"))
redtable$b_D6d <- factor(redtable$b_D6d, levels = c("UP", "NonDE", "DN"))
redtable$b_D1m <- factor(redtable$b_D1m, levels = c("UP", "NonDE", "DN"))
redtable$b_NonD <- factor(redtable$b_NonD, levels = c("UP", "NonDE", "DN"))

ggplot(redtable,
       aes(axis1 = b_D3d, axis2 = b_D6d, axis3 = b_D1m, axis4 = b_NonD)) +
  scale_x_discrete(limits = c("D3d", "D6d", "D1m", "NonD")) +
  geom_alluvium(aes(fill = b_D3d)) +
  geom_stratum(alpha = .5) + 
  geom_text(stat = "stratum", label.strata = TRUE) +
  scale_fill_brewer(type = "qual", palette = "Set2" ) +
  theme(legend.position = "none") +
  ggtitle("gene flow during diapause")
#sessionInfo()


