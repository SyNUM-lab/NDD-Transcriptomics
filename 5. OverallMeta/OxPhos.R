# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load data
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/metaData_all.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("5. overallMeta/GSEA/sets_ordered.RData")
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")
load("Data/CleanData/geneInfo.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)

# Oxidative phosphorylation-related genes
oxphos_genes <- unique(geneInfo$Symbol[geneInfo$GeneID %in% GOgenes[["GO:0006119"]]])
oxphos_ids <- unique(geneInfo$GeneID[geneInfo$Symbol %in% oxphos_genes])

# Get logFCs of OxPhos genes
logFC_matrix[is.na(logFC_matrix)] <- 0
logFC_sel <- logFC_matrix[rownames(logFC_matrix) %in% oxphos_ids,]
plotDF <- gather(as.data.frame(logFC_sel))
plotDF$GeneID <- rep(rownames(logFC_sel), ncol(logFC_sel))
plotDF <- left_join(plotDF, unique(geneInfo[,1:2]), by = "GeneID")

# Distinguish between nuclear- and mitochondrial-encoded genes
mit_genes <- geneInfo[geneInfo$ChrAcc == "NC_012920.1",]$Symbol
plotDF$Mito <- ifelse(plotDF$Symbol %in% mit_genes, "Mitochondrial",
                      "Nuclear")
plotDF$key <- factor(plotDF$key, levels = sets_ordered)

# Make plot
p <- ggplot(plotDF) +
  geom_line(data = plotDF[plotDF$Mito == "Mitochondrial",], 
            aes(x = key, y = value, group = Symbol), alpha = 500/sum(plotDF$Mito == "Mitochondrial"))+
  geom_line(data = plotDF[plotDF$Mito == "Nuclear",], 
            aes(x = key, y = value, group = Symbol), alpha = 500/sum(plotDF$Mito == "Nuclear"))+
  facet_grid(rows = vars(Mito), scale = "free", space = "free") +
  coord_cartesian(ylim = c(-5,5)) +
  xlab("Datasets") +
  ylab(expression(log[2]~"FC")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save plot
ggsave(p, file = "5. OverallMeta/OxPhos/OxPhos.png",
       width = 10, height = 4)


# Find which datasets have the highest peaks
plotDF_fil <- plotDF[plotDF$Mito == "Mitochondrial",]
plotDF_fil <- plotDF_fil %>%
  dplyr::group_by(key) %>%
  summarize(medianVal = median(value))
plotDF_fil <- inner_join(plotDF_fil, metaData_all, by = c("key" = "ID"))


