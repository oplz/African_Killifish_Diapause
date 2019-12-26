## for bulk diapause RNA-Seq data files made by Kallisto

library(sleuth)
library(Biobase)

## working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")   # Set it to your current working directory

##### ===========================================================================================
##### data.frame making 
##### ===========================================================================================

## base directory of the results in a variable
TMF_dir <- "read_mapping"

## name of the Library index file
panelkey<-c("ref_PreD_all5groups")

## to obtain a list of sample IDs in the sub-folder "Nfur" where all data files are
sample_id <- dir(file.path(TMF_dir))

## to obtain a list of paths to the kallisto results
kal_dirs <- sapply(sample_id, function(id) file.path(TMF_dir, id))

## to load an auxillary table that describes the experimental design and the relationship 
## between the kallisto directories and the samples:

s2c <- read.table(file.path(paste0("DE/", panelkey, ".txt")), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path=kal_dirs)


##### ===========================================================================================
##### DE analysis
##### ===========================================================================================

## to load the kallisto processed data into the object
so <- sleuth_prep(s2c, ~ condition)

## to estimate parameters for the sleuth response error measurement model
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

## to perform differential analyis. PreD (Pre-diapause is the Intercept) 
so <- sleuth_wt(so, which_beta = '(Intercept)') 
so <- sleuth_wt(so, which_beta = 'conditionD3d')
so <- sleuth_wt(so, which_beta = 'conditionD6d')
so <- sleuth_wt(so, which_beta = 'conditionD1m')
so <- sleuth_wt(so, which_beta = 'conditionNonD')

models(so)
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

## to generate the Shiny webpage that allows for exploratory data analysis
#sleuth_live(so)

## to save the buildup data (huge file!!)
# sleuth_data_pathway <-paste0(sleuth_data_dir, "/", panelkey, "/", panelkey, "_data.RData")
# save(so, file = sleuth_data_pathway)


##### ===========================================================================================
##### Extract Data from individual condition
##### ===========================================================================================

## to set up a list of conditions.  
conditionlist<-c("(Intercept)",
                 "conditionD3d",
                 "conditionD6d",
                 "conditionD1m",
                 "conditionNonD")

for( i in 1:length(conditionlist)){
  cleaner<-sleuth_results(so, conditionlist[i])

  ## to keep only (target_id), (qval), (beta value)  
  cleaned<-cleaner[,c(1,2,4)]
  
  #convert from ln2 to log2 by multiplying log2(2)/ln(2)
  cleaned[,3]<-cleaned[,3]*1.44269504    
  
  file_names<- paste0("DE/output/differential_gene_list_", conditionlist[i], ".csv")
  write.csv(cleaned, file = file_names)
  rm(cleaner, cleaned)
}


## to extract normalized and filtered TPM from sleuth object 
normalized_TPM<-kallisto_table(so, use_filtered=TRUE, normalized = TRUE, include_covariates = TRUE)
write.csv(normalized_TPM, file="DE/output/normalized_TPM.csv")

rm(list=setdiff(ls(), "panelkey"))


##### ===========================================================================================
##### post-analysis list making
##### ===========================================================================================

#### (1) to make a list of raw read counts (TPM) from all libraries

## to use Lib01 (PreD_01) as the template 
seedTPMlist<-read.csv(paste0("read_mapping/01/abundance.tsv"), 
                header= TRUE, colClasses = c("character", "NULL","NULL", "NULL", "numeric"),
                sep = "\t")

## to add other Libs by order (all with the same gene order in the reference transcriptome; no need to reorder)
TPM_Lib_list<-c("02", "03",
                "16", "15", "14",
                "13", "12", "11",
                "10", "09", "08",
                "04", "05", "06")

for (i in 1:length(TPM_Lib_list)) {
  newaddTPMlist<-read.csv(paste0("read_mapping/", TPM_Lib_list[i], "/abundance.tsv"), 
                header= TRUE, colClasses = c("NULL", "NULL","NULL", "NULL", "numeric"),
                sep = "\t")
  seedTPMlist<-cbind(seedTPMlist, newaddTPMlist)
  rm(newaddTPMlist)
}

rownames(seedTPMlist)<-seedTPMlist$target_id
seedTPMlist$target_id<-NULL
colnames(seedTPMlist)<-c("TPM_PreD_01", "TPM_PreD_02", "TPM_PreD_03",
                      "TPM_D3d_11", "TPM_D3d_12", "TPM_D3d_24",
                      "TPM_D6d_14", "TPM_D6d_15", "TPM_D6d_16",
                      "TPM_D1m_18", "TPM_D1m_19", "TPM_D1m_26",
                      "TPM_NonD_21", "TPM_NonD_22", "TPM_NonD_23")

seedTPMlist<-seedTPMlist[order(rownames(seedTPMlist)),]  ## reorder the genes in alphabetrical order 

write.csv(seedTPMlist, file="DE/output/RAW_TPM_combined.csv")
rm(TPM_Lib_list,i)

#### (2) to make a list of DE information in all libraries

## to load PreD(Intercept) as the template, but only keep the Gene ID (the reference for DE).  
seedDElist<-read.csv(paste0("DE/output/differential_gene_list_(Intercept).csv"), 
                     header= TRUE, colClasses = c("NULL", "character", "NULL", "NULL"),
                     sep = ",")

rownames(seedDElist)<-seedDElist$target_id
seedDElist$target_id<-NULL
seedDElist<-seedDElist[order(rownames(seedDElist)),]   ## reorder the genes in alphabetrical order 

## to append other conditions by order
DE_condition_list<-c("D3d", "D6d","D1m", "NonD")

for (i in 1:length(DE_condition_list)) {
  newaddDElist<-read.csv(paste0("DE/output/differential_gene_list_condition", DE_condition_list[i], ".csv"), 
                      header= TRUE, colClasses = c("NULL", "character", "numeric", "numeric"),
                      sep = ",")
  rownames(newaddDElist)<-newaddDElist$target_id
  newaddDElist$target_id<-NULL
  newaddDElist<-newaddDElist[order(rownames(newaddDElist)),]   ## reorder the genes in alphabetrical order  
  seedDElist<-cbind(seedDElist, newaddDElist)
  rm(newaddDElist)
}

colnames(seedDElist)<-c("qval_D3d", "b_D3d", "qval_D6d", "b_D6d", "qval_D1m", "b_D1m", "qval_NonD", "b_NonD")
write.csv(seedDElist, file="DE/output/DE_combined.csv")
rm(DE_condition_list, i)

#### (3) to combine both DE and TPM lists into one master list
combined<-cbind(seedDElist, seedTPMlist)
write.csv(combined, file="DE/output/Master_combined.csv")
rm(list=ls())

