# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Perform GSEA for each study
for (i in 1:length(topList)){

  # Maked sorted list based on logFC
  gsea_input <- sort(setNames(topList[[i]]$log2FC,topList[[i]]$GeneID), 
                     decreasing = TRUE)
  
  # Perform GSEA
  set.seed(123)
  GOtest <- gseGO(
    geneList = gsea_input,
    keyType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    ont =  "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500)
  
  # Get results
  GOresults <- GOtest@result
  
  # Collect P values and NES
  if (i == 1){
    pvalues <- GOresults[,c("ID", "Description", "pvalue")]
    colnames(pvalues) <- c("ID", "Description", names(topList)[i])
    NES <- GOresults[,c("ID", "Description", "NES")]
    colnames(NES) <- c("ID", "Description", names(topList)[i])
  } else{
    temp_p <- GOresults[,c("ID", "pvalue")]
    colnames(temp_p) <- c("ID", names(topList)[i])
    pvalues <- full_join(pvalues, temp_p, by = c("ID" = "ID"))
    
    temp_n <- GOresults[,c("ID", "NES")]
    colnames(temp_n) <- c("ID", names(topList)[i])
    NES <- full_join(NES, temp_n, by = c("ID" = "ID"))
  }
  
}

save(pvalues, NES, file = "4. GSEA/GSEA_GO/GSEAresults_GO.RData")



