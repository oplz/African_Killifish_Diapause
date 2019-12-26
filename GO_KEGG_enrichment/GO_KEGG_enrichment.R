## for pathway enrichment of GO and KEGG 

## working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory



##### ===========================================================================================
##### Parameters
##### ===========================================================================================

## location of the input gene list 
# foldername = "BLE_CBX7_N+C"
# foldername = "DE_heatmap"
# foldername = "DE_p001"
# foldername = "PCA_components"
# foldername = "DEG_cluster6"
foldername = "working"

## GO term cutoff value.  Use 1 to have no cutoff (for GO term heatmap)
GO_P_cutoff_value = 0.01 
KEGG_P_cutoff_value = 0.05

## Minimum number of genes for a term to filter
GO_mingenes = 5
KEGG_mingenes = 2

## Relative enrichment filter
GO_relenrich = 1
KEGG_relenrich = 2

## ontology MF, BP, CC
ontolg = "BP"


##### ===========================================================================================
##### GO term richment
##### ===========================================================================================
library("GOstats")

## to check if the saving directory exists, if not, build one, using "foldername".
if (file.exists(paste0("GO_KEGG_enrichment/output/", foldername))){
  warningmsg<-paste0("Data outputing folder => [GO_KEGG_enrichment/output/", foldername, "] already exists!")
  print(warningmsg)
} else {
  dir.create(paste0("GO_KEGG_enrichment/output/", foldername))
}

## to read all lists except rhe background (BG.txt)
pairlist<-list.files(paste0("gene_lists/", foldername), pattern = ".txt", recursive = TRUE)
pairlist<-pairlist[pairlist != "BG.txt"]

## All GO terms
frame = read.table(file ="GO_KEGG_enrichment/reference/GO_killifish-human_best-hits.txt", header = T, colClasses=c(rep("factor",3)))
## List of universe genes. Background.
universe = read.table(file = paste0("gene_lists/", foldername, "/BG.txt"), header = T)
mingenes = GO_mingenes
relenrich = GO_relenrich

## to put data into a GOFrame object
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
goFrame=GOFrame(goframeData,organism="Human")

## to cast this object to a GOAllFrame object will tap into the GO.db package and populate this object with the implicated GO2All mappings for you
goAllFrame=GOAllFrame(goFrame)

## generate geneSetCollection objects
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

## to process the universe list
universe = universe$id
universe = lapply(universe, as.character)
universe = unlist(universe)

for (cycling in 1:length(pairlist)) {
## to process the gene list of interest
  genes = read.table(file=paste0("gene_lists/", foldername, "/", pairlist[cycling]), header = T)
  genes = genes$id
  genes = lapply(genes, as.character)
  genes = unlist(genes)

  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params", 
                               geneSetCollection=gsc, geneIds = genes, 
                               universeGeneIds = universe,
                               ontology = ontolg,
                               pvalueCutoff = GO_P_cutoff_value,
                               conditional = FALSE, 
                               testDirection = "over") # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over"

## to call hyperGTest to do the test
  Over <- hyperGTest(params)

  ## to calculate enrichment and add it to data frame.
  ## Relative enrichment factor (E-value) for a GO term = (count/size)/(size/universe)  
  enrichment = (summary(Over)[5]$Count / summary(Over)[6]$Size) / (summary(Over)[6]$Size / length(universe))
  
  ## to create a new frame
  SummaryOver = data.frame(summary(Over), enrichment)

  ## to filter the Over variable on parameters other than P-value
  ## to filter the summary of OVER with size of the term, at least 2 genes for a go term
  FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$enrichment >= relenrich),]

  ## to adjust p value for multile correction
  padj = p.adjust(FilteredSummaryOver$Pvalue, "BH")
  
  ## to add padj to the data frame
  FinalSummaryOver = data.frame(FilteredSummaryOver, padj)
  
  ## to write to a file
  write.table(FinalSummaryOver, paste0("GO_KEGG_enrichment/output/", foldername, "/", "GO_", pairlist[cycling]), quote = F, row.names = F, sep = "\t")
  
## to isolate indexes for the go terms in final results
  ind.GO <- is.element(names(Over@goDag@nodeData@data), FinalSummaryOver$GOBPID) 
  selected.GO <- Over@goDag@nodeData@data[which(ind.GO)]

## to get a go terms and genes in a new vaiable for all the terms in the results of enrichment
  goTerms <- lapply(selected.GO, function(x) x$geneIds)
  names(goTerms) <- names(Over@goDag@nodeData@data)[ind.GO]

## to create the files containing the genes in individual GO terms
  for (i in 1:length(goTerms)){
    test = as.data.frame(do.call(rbind, goTerms[i]))
    write.table(test, paste0("GO_KEGG_enrichment/output/", foldername, "/", "GO_genes_", pairlist[cycling]), quote = F, col.names = F, append = T) 
    rm(test)
  }
}


## to create a data list to keep and erase the rest
safelist<-c("foldername", "pairlist", "KEGG_relenrich", "KEGG_mingenes", "KEGG_P_cutoff_value")
rm(list=setdiff(ls(), safelist))
   

##### ===========================================================================================
##### KEGG enrichment
##### ===========================================================================================

## all KEGG terms
frame = read.table(file ="GO_KEGG_enrichment/reference/KEGG_killifish-human_best-hits.txt", header = T, sep = "\t", colClasses=c(rep("factor",2)))

## to obtain information of kegg_id and gene_id
keggframeData = data.frame(frame$kegg_id, frame$gene_id)

## list of universe genes. background.
universe = read.table(file = paste0("gene_lists/", foldername, "/BG.txt"), header = T)
mingenes = KEGG_mingenes
relenrich = KEGG_relenrich

library("KEGG.db")
## to put data into a GOFrame object
keggFrame=KEGGFrame(keggframeData, organism="Human")

## to generate geneSetCollection objects
gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())

## to process the universe list universe
universe = universe$id
universe = lapply(universe, as.character)
universe = unlist(universe)


## to process the list of interest
for (cycling in 1:length(pairlist)) {
  ## to process the gene list of interest
  genes = read.table(file=paste0("gene_lists/", foldername, "/", pairlist[cycling]), header = T)
  genes = genes$id
  genes = lapply(genes, as.character)
  genes = unlist(genes)
  
  kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based KEGG enrichment analysis", 
                                  geneSetCollection=gsc, 
                                  geneIds = genes, 
                                  universeGeneIds = universe,
                                  pvalueCutoff = KEGG_P_cutoff_value,
                                  testDirection = "over") # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over"
  
  ## to call hyperGTest to do the test
  kOver <- hyperGTest(kparams)
  #head(summary(kOver))
  
  ## to adjust p value for multile correction
  kpadj = p.adjust(summary(kOver)[2]$Pvalue, "BH")
  
  ## to calculate enrichment and add it to data frame.
  #Enrichmen = (count/size)/(size/13637) where 13637 is all the genes considered for selection anakysis
  enrichment = (summary(kOver)[5]$Count / summary(kOver)[6]$Size) / (summary(kOver)[6]$Size / length(universe))
  
  
  ## to create a new frame
  SummarykOver = data.frame(summary(kOver), enrichment)
  
  ## to filter the Over variable on parameters other than P-value
  ## to filter the summary of OVER with size of the term, at least 2 genes for a go term
  FilteredSummarykOver = SummarykOver[(SummarykOver$Count >= mingenes & SummarykOver$enrichment >= relenrich),]
  
  # to adjust p value for multile correction
  kpadj = p.adjust(FilteredSummarykOver$Pvalue, "BH")
  
  # to add padj to the data frame
  FinalSummarykOver = data.frame(FilteredSummarykOver, kpadj)
  
  # to write to a file
  write.table(FinalSummarykOver, paste0("GO_KEGG_enrichment/output/", foldername, "/", "KEGG_", pairlist[cycling], ".txt"), quote = F, row.names = F, sep = "\t")
  
}

rm(list = ls())
