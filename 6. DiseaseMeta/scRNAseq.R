# Data available from:
# https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(Seurat)
library(sparseMatrixStats)
library(patchwork)

selGenes <- c("ENSG00000132470", "ENSG00000070087",
              "ENSG00000165512", "ENSG00000068784")
selNames <- c("ITGB4", "PFN2", "ZNF22", "SRBD1")


files <- list.files("Data/BrainAtlas")
mRes <- matrix(NA, nrow = length(files), ncol = length(selGenes))
nCel <- matrix(NA, nrow = length(files), ncol = length(selGenes))
for (f in 1:length(files)){
  # Read data
  exprData <- read_rds(paste0("Data/BrainAtlas/", files[f]))
  
  for (i in 1:length(selGenes)){
    # Get expression
    expr <- exprData@assays$RNA@data[selGenes[i],]
    
    # Get median
    nCel[f,i] <- sum(expr != 0)/length(expr)
    mRes[f,i] <- mean(expr[expr!=0])
  }

}
rownames(mRes) <- str_remove(files, ".rds")
rownames(nCel) <- str_remove(files, ".rds")
colnames(mRes) <- selNames
colnames(nCel) <- selNames

# Save data
save(mRes, nCel, file = "6. DiseaseMeta/scRNAseq/markerResults.RData")

# Load data
load("6. DiseaseMeta/scRNAseq/markerResults.RData")
mRes_scaled <- apply(mRes,2,function(x) x/max(x))

# Format data for plotting
plotDF <- gather(as.data.frame(mRes_scaled))
plotDF$cluster <- rep(rownames(mRes_scaled), ncol(mRes_scaled))
plotDF$nCel <- gather(as.data.frame(nCel))$value
plotDF$cluster[plotDF$cluster == "Deep-layer corticothalamic and 6b"] <- "Deep-layer corticothalamic "


nonneuron <- c("Astrocyte", "Bergman glia", "Choroid plexus",
               "Committed oligendrocyte precursor", "Ependymal",
               "Fibroblast", "Microglia", "Oligodendrocyte precursor", "Oligodendrocyte",
               "Vascular")

plotDF$Group <- ifelse(plotDF$cluster %in% nonneuron, "Non-neuronal", "Neuronal")
plotDF$key <- factor(plotDF$key, levels = rev(c("PFN2", "ZNF22", "SRBD1", "ITGB4")))


# Make plot
p <- ggplot(plotDF) +
  geom_point(aes(x = cluster, y = key, size = nCel, color = value)) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  scale_color_gradient(low = "#D9D9D9", high = "#CB181D")+
  #scale_color_gradient(low="#7570B3", high="#E7298A")+
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal() +
  labs(size = "Cell\nproportion", color = "Scaled mean\nexpression") +
  scale_x_discrete(labels = scales::label_wrap(18)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 12))

# Save plot
ggsave(p, file = "6. DiseaseMeta/scRNAseq/cellExpr2.png", width = 10, height = 3)


# Make legend plot
legendPlot <- ggplot(plotDF) +
  geom_point(aes(x = cluster, y = key, size = nCel, color = value)) +
  scale_color_gradient(low = "#D9D9D9", high = "#CB181D")+
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  labs(size = "Cell\nproportion", color = "Scaled mean\nexpression") +
  scale_x_discrete(labels = scales::label_wrap(18)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "6. DiseaseMeta/scRNAseq/legend1.png", width = 8, height = 8)


# Export figure source data
sourceData <- plotDF[,c("key", "cluster", "Group", "value", "nCel")]
colnames(sourceData) <- c("Gene", "CellType", "Cluster", "Expr", "Proportion")
write.csv(sourceData, file = "6. DiseaseMeta/SourceData_Figure3C.csv",
          row.names = FALSE, quote = FALSE)
