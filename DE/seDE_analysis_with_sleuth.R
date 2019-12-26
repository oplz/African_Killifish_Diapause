## for single embryo diapause RNA-Seq data files made by Kallisto

library(sleuth)
library(Biobase)

## to set the "Diapause_transcritome" as the working directory
setwd("/Users/Chi-Kuo/Dropbox/Diapause_transcriptome")    # Set it to your current working directory

##### ===========================================================================================
##### data.frame making 
##### ===========================================================================================

## base directory of the results in a variable
TMF_dir <- "se_read_mapping"

## the name of the Library index file
panelkey<-c("ref_BLE")

## to obtain a list of sample IDs in the sub-folder "Nfur" where all data files are
sample_id <- dir(file.path(TMF_dir))

## to obtain a list of paths to the kallisto results indexed by the sample IDs is collated with
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
so <- sleuth_wt(so, which_beta = 'conditionBLEBB')

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
conditionlist<-c("(Intercept)", "conditionBLEBB")

for( i in 1:length(conditionlist)){
  cleaner<-sleuth_results(so, conditionlist[i])

  ## to keep only column 1 (target_id), 3 (qval), and 4 (beta value)  
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
##### post-analysis lists making
##### ===========================================================================================

#### (1) to make a list of raw read counts (TPM) from all libraries

## to use Lib01 (PreD_01) as the template 
seedTPMlist<-read.csv(paste0("se_read_mapping/0D239/abundance.tsv"), 
                header= TRUE, colClasses = c("character", "NULL","NULL", "NULL", "numeric"),
                sep = "\t")

## to add other Libs by order (all with the same gene order in the reference transcriptome; no need to reorder)

TPM_Lib_list<-c("0D2310", "0D2317", "0D2318", "0D2325", "0D2326",
                "0D2333", "0D2334", "0D2335", "0D2336", "0D2337",
                "BLEBB5", "BLEBB6", "BLEBB7", "BLEBB8", "BLEBB13",
                "BLEBB14", "BLEBB16", "BLEBB21", "BLEBB22", "BLEBB23",
                "BLEBB24", "BLEBB29", "BLEBB30", "BLEBB31", "BLEBB32",
                "BLEBB41", "BLEBB42"
)

for (i in 1:length(TPM_Lib_list)) {
  newaddTPMlist<-read.csv(paste0("se_read_mapping/", TPM_Lib_list[i], "/abundance.tsv"), 
                header= TRUE, colClasses = c("NULL", "NULL","NULL", "NULL", "numeric"),
                sep = "\t")
  seedTPMlist<-cbind(seedTPMlist, newaddTPMlist)
  rm(newaddTPMlist)
}

rownames(seedTPMlist)<-seedTPMlist$target_id
seedTPMlist$target_id<-NULL
colnames(seedTPMlist)<-c("0D239","0D2310", "0D2317", "0D2318", "0D2325", "0D2326",
                         "0D2333", "0D2334", "0D2335", "0D2336", "0D2337",
                         "BLEBB5", "BLEBB6", "BLEBB7", "BLEBB8", "BLEBB13",
                         "BLEBB14", "BLEBB16", "BLEBB21", "BLEBB22", "BLEBB23",
                         "BLEBB24", "BLEBB29", "BLEBB30", "BLEBB31", "BLEBB32",
                         "BLEBB41", "BLEBB42"
                         )

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
DE_condition_list<-c("BLEBB")

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
colnames(seedDElist)<-c("qval_BLEBB", "b_BLEBB")
write.csv(seedDElist, file="DE/output/DE_combined.csv")
rm(DE_condition_list, i)

#### (3) to combine both DE and TPM lists into one master list
combined<-cbind(seedDElist, seedTPMlist)
write.csv(combined, file="DE/output/Master_combined.csv")
rm(list=ls())

