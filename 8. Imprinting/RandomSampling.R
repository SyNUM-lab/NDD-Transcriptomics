# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")
load("Data/CleanData/geneInfo.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/topList.RData")

# Load packages
library(tidyverse)
library(readxl)

# Get upregulated, downregulated, and differentially expressed genes
up_genes <- list()
down_genes <- list()
both_genes <- list()
all_genes <- list()
for (t in 1:length(topList)){
  testList <- topList[[t]][!is.na(topList[[t]]$`p-value`),]
  up_genes[[t]] <- testList$GeneID[(testList$`p-value` < 0.05) & (testList$log2FC > 0)]
  down_genes[[t]] <- testList$GeneID[(testList$`p-value` < 0.05) & (testList$log2FC < 0)]
  both_genes[[t]] <- testList$GeneID[(testList$`p-value` < 0.05)]
  all_genes[[t]] <- testList$GeneID
}

# Extract mean expression
for (i in 1:length(topList)){
  if (i > 1) {
    meanExpr_temp <- topList[[i]][,c(1,2)]
    meanExpr <- full_join(meanExpr, meanExpr_temp, by = c("GeneID" = "GeneID"))
  } else{
    meanExpr <- topList[[i]][,c(1,2)]
  }
}
meanExpr$GeneID <- as.character(meanExpr$GeneID)
colnames(meanExpr) <- c("GeneID", names(topList))
meanExpr_matrix <- as.matrix(meanExpr[,-1])
rownames(meanExpr_matrix) <- meanExpr$GeneID

# Get quantile of median expression
meanExpr_matrix[is.na(meanExpr_matrix)] <- 0
medianExpr <- rev(sort(apply(meanExpr_matrix,1,median)))
Q1 <- names(medianExpr[1:round(length(medianExpr)*0.25)])
Q2 <- names(medianExpr[(round(length(medianExpr)*0.25)+1):round(length(medianExpr)*0.50)])
Q3 <- names(medianExpr[(round(length(medianExpr)*0.50)+1):round(length(medianExpr)*0.75)])
Q4 <- names(medianExpr[(round(length(medianExpr)*0.75)+1):length(medianExpr)])

# Get expression profile of imprinted genes
imprintDF <- read_xlsx("Data/ImprintedGenes.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]
imprinted_genes1 <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
imprinted_genes1 <- imprinted_genes1[!is.na(imprinted_genes1)]
imprinted_genes1 <- unique(as.character(geneInfo$GeneID[geneInfo$Symbol %in% imprinted_genes1]))

# Get number of imprinted genes per expression quantile
Q1imp <- sum(which(names(medianExpr) %in% imprinted_genes1)/length(medianExpr) <= 0.25)
Q2imp <- sum(which(names(medianExpr) %in% imprinted_genes1)/length(medianExpr) <= 0.50) - Q1imp
Q3imp <- sum(which(names(medianExpr) %in% imprinted_genes1)/length(medianExpr) <= 0.75) - Q1imp - Q2imp
Q4imp <- sum(which(names(medianExpr) %in% imprinted_genes1)/length(medianExpr) <= 1) - Q1imp - Q2imp - Q3imp


# Perform permuations analysis
nperm <- 1000
results_up <- matrix(NA,nrow = length(both_genes), ncol = nperm)
results_down <- matrix(NA,nrow = length(both_genes), ncol = nperm)
results_both <- matrix(NA,nrow = length(both_genes), ncol = nperm)
set.seed(123)
for (p in 1:nperm){
  sel_genes <- c(Q1[sample(1:length(Q1),Q1imp)],
                 Q2[sample(1:length(Q2),Q2imp)],
                 Q3[sample(1:length(Q3),Q3imp)],
                 Q4[sample(1:length(Q4),Q4imp)])
  
  for (d in 1:length(both_genes)){
    
    # get all genes
    all_genes1 <- all_genes[[d]]
    
    # imprinted genes
    imprinted_genes <- intersect(all_genes1, sel_genes)
    
    # get non-imprinted genes
    nimprinted_genes <- setdiff(all_genes1, imprinted_genes)
    nimprinted_genes <- nimprinted_genes[!is.na(nimprinted_genes)]
    
    # Evaluate upregulated genes:
    
    # get upregulated genes
    up_genes1 <- up_genes[[d]]
    
    # get non-upregulated genes
    nup_genes1 <- setdiff(all_genes1, up_genes1)
    
    
    up_imp <- intersect(up_genes1,imprinted_genes)
    nup_imp <- intersect(nup_genes1,imprinted_genes)
    up_nimp <- intersect(up_genes1, nimprinted_genes)
    nup_nimp <- intersect(nup_genes1, nimprinted_genes)
    
    m <- matrix(c(length(up_imp), 
                  length(nup_imp), 
                  length(up_nimp),
                  length(nup_nimp)),nrow = 2)
    
    output <- fisher.test(m)
    results_up[d,p] <- output$estimate
    
    # Evaluate downregulated genes:
    
    # get downpregulated genes
    down_genes1 <- down_genes[[d]]
    
    # get non-downregulated genes
    ndown_genes1 <- setdiff(all_genes1, down_genes1)
    
    
    down_imp <- intersect(down_genes1,imprinted_genes)
    ndown_imp <- intersect(ndown_genes1,imprinted_genes)
    down_nimp <- intersect(down_genes1, nimprinted_genes)
    ndown_nimp <- intersect(ndown_genes1, nimprinted_genes)
    
    m <- matrix(c(length(down_imp), 
                  length(ndown_imp), 
                  length(down_nimp),
                  length(ndown_nimp)),nrow = 2)
    
    output <- fisher.test(m)
    results_down[d,p] <- output$estimate
    
    
    # Evaluate up- and downregulated genes:
    
    # get DEd genes
    both_genes1 <- both_genes[[d]]
    
    # get non-DEd genes
    nboth_genes1 <- setdiff(all_genes1, both_genes1)
    
    
    both_imp <- intersect(both_genes1,imprinted_genes)
    nboth_imp <- intersect(nboth_genes1,imprinted_genes)
    both_nimp <- intersect(both_genes1, nimprinted_genes)
    nboth_nimp <- intersect(nboth_genes1, nimprinted_genes)
    
    m <- matrix(c(length(both_imp), 
                  length(nboth_imp), 
                  length(both_nimp),
                  length(nboth_nimp)),nrow = 2)
    
    output <- fisher.test(m)
    results_both[d,p] <- output$estimate
  }
  
  
}

# Format data
for (c in 1:ncol(results_both)){
  results_both[,c] <- sort(results_both[,c])
  results_up[,c] <- sort(results_up[,c])
  results_down[,c] <- sort(results_down[,c])
}

# Save results
save(results_both, results_down, results_up,
     file = "8. Imprinting/PermResults1.RData")

