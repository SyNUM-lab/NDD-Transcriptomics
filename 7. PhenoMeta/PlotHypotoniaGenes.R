# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(clusterProfiler)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/geneInfo.RData")

# Select gene and phenotype of interest
selGenes <- c("64968", "51227", "2114", "7327")
selPheno <- "Hypotonia"
selSamples <- metaData_all$ID[metaData_all[,selPheno] == 1]

# Prepare data for plotting
plotDF_all <- NULL
for (i in 1:length(selGenes)){
  selGene <- selGenes[i]
  name <- geneInfo$Symbol[geneInfo$GeneID == selGene]
  
  # Format logFC and SE matrices
  logFC_matrix1 <- logFC_matrix
  logFC_matrix1[is.na(logFC_matrix1)] <- 0
  SE_matrix1 <- SE_matrix
  SE_matrix1[is.na(SE_matrix1)] <- 0
  
  # Prepare data for plotting
  plotDF <- data.frame(logFC = logFC_matrix1[which(rownames(logFC_matrix1) == selGene),],
                       Lower = logFC_matrix1[which(rownames(logFC_matrix1) == selGene),] - 1.96*SE_matrix1[which(rownames(logFC_matrix1) == selGene),],
                       Upper = logFC_matrix1[which(rownames(logFC_matrix1) == selGene),] + 1.96*SE_matrix1[which(rownames(logFC_matrix1) == selGene),],
                       Pvalue = pvalue_matrix[which(rownames(pvalue_matrix) == selGene),],
                       Pheno = factor(ifelse(colnames(logFC_matrix1) %in% selSamples, selPheno, paste0("No ", selPheno)),
                                      levels = c(selPheno, paste0("No ", selPheno))),
                       Dataset = colnames(logFC_matrix1),
                       GeneID = selGene,
                       GeneName = name)
  plotDF$Sig <- ifelse((plotDF$Lower > 0) | (plotDF$Upper < 0), "Yes", "No")
  
  
  # Combine with meta data
  plotDF <- inner_join(plotDF, metaData_all, by = c("Dataset" = "ID"))
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

# Make plot
p <- ggplot(plotDF_all) +
  geom_point(aes(x = Dataset, y = logFC, 
                 color = logFC)) +
  geom_segment(aes(x = Dataset, xend = Dataset, y = Lower, yend = Upper, 
                   color = logFC)) +
  facet_grid(cols = vars(Pheno), rows = vars(GeneName),
             scale = "free", space = "free_x") +
  scale_color_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                        trans = "pseudo_log") +
  ylab(expression(log[2]~"FC")) +
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  xlab("Datasets") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, hjust = 1),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        strip.text.x = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 12))
# Save plot
ggsave(p, file = "7. PhenoMeta/Gene_Pheno/Hypotonia.png", width = 8, height = 3)
