
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)

# Read and prepare data
datasets <- list.files("Data/Datasets")
topList <- list()

for (d in 1:length(datasets)){
  topfiles <- list.files(paste0("Data/Datasets/", datasets[d]), pattern = "topTable")
  
  if (length(topfiles) >0 ){
    
    if (length(topfiles) > 1){
      
      for (t in 1:length(topfiles)){
        name <- paste0(datasets[d],"-",str_remove(str_remove(topfiles[t], "topTable"), "_.*"))
        topList[[name]] <- data.table::fread(paste0("Datasets/", datasets[d],"/",topfiles[t]))
      }
      
    } else{
      name <- datasets[d]
      topList[[name]] <- data.table::fread(paste0("Datasets/", datasets[d],"/",topfiles))
    }
  }
}


# Check if information is missing from topList or MetaData
meta <- readxl::read_excel("Data/CleanData/MetaData_clean.xlsx")
setdiff(meta$ID, names(topList))
setdiff(names(topList), meta$ID)

# Save topList
save(topList, file = "Data/CleanData/topList.RData")

