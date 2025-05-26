# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(readxl)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Load data
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/geneInfo.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)

# Get upper and low C.I.
logFC_matrix[is.na(logFC_matrix)] <- 0
SE_matrix[is.na(SE_matrix)] <- 10
upper <- logFC_matrix + 1.96*SE_matrix
lower <- logFC_matrix - 1.96*SE_matrix

# Get number of significant datasets
testDF <- data.frame(GeneID = rownames(logFC_matrix),
                     value = rowSums((lower > 0) | (upper < 0)))
testDF <- inner_join(testDF, unique(geneInfo[,1:2]), by = c("GeneID" = "GeneID"))

# Get imprinted genes: https://www.geneimprint.com/site/genes-by-species
imprintDF <- read_xlsx("Data/ImprintedGenes.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]
imprinted_genes1 <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
imprinted_genes1 <- imprinted_genes1[!is.na(imprinted_genes1)]

# Select imprinted genes only
testDF_imp <- testDF[testDF$Symbol %in% imprinted_genes1,]


# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(readxl)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Load data
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/geneInfo.RData")
load("5. overallMeta/GSEA/sets_ordered.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)

# Format logFC and SE matrices
logFC_matrix1 <- logFC_matrix
logFC_matrix1[is.na(logFC_matrix1)] <- 0
SE_matrix1 <- SE_matrix
SE_matrix1[is.na(SE_matrix1)] <- 0

# Select gene
selName <- "MEST"
selGene <- unique(geneInfo$GeneID[geneInfo$Symbol == selName])

# Prepare data for plotting
plotDF <- data.frame(logFC = logFC_matrix1[which(rownames(logFC_matrix1) == selGene),],
                     Lower = logFC_matrix1[which(rownames(logFC_matrix1) == selGene),] - 1.96*SE_matrix1[which(rownames(logFC_matrix1) == selGene),],
                     Upper = logFC_matrix1[which(rownames(logFC_matrix1) == selGene),] + 1.96*SE_matrix1[which(rownames(logFC_matrix1) == selGene),],
                     Pvalue = pvalue_matrix[which(rownames(pvalue_matrix) == selGene),],
                     Dataset = colnames(logFC_matrix1))

plotDF$Sig <- ifelse((plotDF$Lower > 0) | (plotDF$Upper < 0) , "Yes", "No")
plotDF <- inner_join(plotDF, metaData_all, by = c("Dataset" = "ID"))
plotDF$Dataset <- factor(plotDF$Dataset,
                         levels = sets_ordered)

# Make plot
colors <- c("#BDBDBD", "#6A51A3")
p <- ggplot(plotDF) +
  geom_point(aes(x = Dataset, y = logFC, color = Sig), size = 1) +
  geom_segment(aes(x = Dataset, xend = Dataset, y = Lower, yend = Upper, color = Sig)) +
  ylab(expression(log[2]~"FC")) +
  xlab("Datasets") +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ggtitle(selName) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        strip.text = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 12))

# Save
ggsave(p, file = paste0("8. Imprinting/Figures/",selName, "_logFCs.png"), width = 9, height = 2)
