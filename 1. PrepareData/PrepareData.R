# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/Analysis")

# Load packages
library(tidyverse)
library(biomaRt)

# Load data
load("Data/CleanData/topList.RData")
metaData <- readxl::read_excel("Data/CleanData/MetaData_clean.xlsx")
geneInfo <- data.table::fread("Data/CleanData/Human.GRCh38.p13.annot.tsv")

################################################################################

# Format gene information

################################################################################

# Format chromosome name
chrAnn <- data.frame(ChrAcc = c("NC_000001.11", "NC_000002.12", "NC_000003.12",
                                "NC_000004.12", "NC_000005.10", "NC_000006.12",
                                "NC_000007.14", "NC_000008.11", "NC_000009.12",
                                "NC_000010.11", "NC_000011.10", "NC_000012.12",
                                "NC_000013.11", "NC_000014.9", "NC_000015.10",
                                "NC_000016.10", "NC_000017.11", "NC_000018.10",
                                "NC_000019.10", "NC_000020.11", "NC_000021.9",
                                "NC_000022.11", "NC_000023.11", "NC_000024.10"),
                     ChrName = c("1", "2", "3", 
                                 "4", "5", "6",
                                 "7", "8", "9",
                                 "10", "11", "12",
                                 "13", "14", "15",
                                 "16", "17", "18",
                                 "19", "20", "21",
                                 "22","X", "Y"))

geneInfo <- left_join(geneInfo, chrAnn, by = c("ChrAcc" = "ChrAcc"))
geneInfo <- unique(geneInfo[,c("GeneID", "Symbol", "Description", "Synonyms", "GeneType",
                        "EnsemblGeneID", "ChrName", "ChrAcc", "ChrStart", "ChrStop",
                        "Orientation", "Length")])

# Save gene information
save(geneInfo, file = "Data/CleanData/geneInfo.RData")


################################################################################

# Format statistics

################################################################################

# Extract P values and logFCs
for (i in 1:length(topList)){
  if (i > 1) {
    logFCs_temp <- topList[[i]][,c(1,3)]
    logFCs <- full_join(logFCs, logFCs_temp, by = c("GeneID" = "GeneID"))
    
    SEs_temp <- topList[[i]][,c(1,4)]
    SEs <- full_join(SEs, SEs_temp, by = c("GeneID" = "GeneID"))
    
    pvalues_temp <- topList[[i]][,c(1,5)]
    pvalues <- full_join(pvalues, pvalues_temp, by = c("GeneID" = "GeneID"))
  } else{
    logFCs <- topList[[i]][,c(1,3)]
    SEs <- topList[[i]][,c(1,4)]
    pvalues <- topList[[i]][,c(1,5)]
  }
}

# Format logFCs
logFCs$GeneID <- as.character(logFCs$GeneID)
colnames(logFCs) <- c("GeneID", names(topList))
logFC_matrix <- as.matrix(logFCs[,-1])
rownames(logFC_matrix) <- logFCs$GeneID

# Format SEs
SEs$GeneID <- as.character(SEs$GeneID)
colnames(SEs) <- c("GeneID", names(topList))
SE_matrix <- as.matrix(SEs[,-1])
rownames(SE_matrix) <- SEs$GeneID

# Format P values
pvalues$GeneID <- as.character(pvalues$GeneID)
colnames(pvalues) <- c("GeneID", names(topList))
pvalue_matrix <- as.matrix(pvalues[,-1])
rownames(pvalue_matrix) <- pvalues$GeneID

# Save statistics
save(logFC_matrix, SE_matrix, pvalue_matrix, file = "Data/CleanData/statistics_matrix.RData")



