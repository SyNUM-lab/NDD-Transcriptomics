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
library(rrvgo)
library(RCy3)
library(igraph)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")
load("5. OverallMeta/GSEA/reducedTerms.RData")

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Jaccard Index
JI <- function(x,y){length(intersect(x,y))/length(union(x,y))}

# Prepare data
selTerms <- reducedTerms[,c("parent", "parentTerm", "score")]
selTerms <- selTerms %>%
  group_by(parent) %>%
  mutate(score = max(score))
rownames(selTerms) <- NULL
selTerms <- unique(selTerms)
plotTerms <- head(selTerms$parent,30)
plotNames <- firstup(head(selTerms$parentTerm,30))
plotNames[nchar(plotNames)>50] <- paste0(substring(plotNames[nchar(plotNames)>50],1,47),"...")

# Make a matrix that shows pairwise Jaccard Index
GOgenes_fil <- GOgenes[plotTerms]
graph_matrix <- matrix(NA, nrow = length(plotTerms), ncol = length(plotTerms))
for (i in 1:length(plotTerms)){
  graph_matrix[i,] <- unlist(lapply(GOgenes_fil,function(x){JI(x,GOgenes_fil[[i]])}))
}
rownames(graph_matrix) <- plotNames
colnames(graph_matrix) <- plotNames

# make a graph from this matrix
graph_matrix_fil <- graph_matrix
g <- igraph::graph_from_adjacency_matrix(graph_matrix_fil, 
                                         mode = "lower", 
                                         weighted = "Jaccard Index", 
                                         diag = FALSE)
# Export network to cytoscape
createNetworkFromIgraph(
  g,
  title = "GO terms",
  collection = "Overall meta-analysis",
)

