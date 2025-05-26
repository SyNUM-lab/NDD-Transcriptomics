
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(clusterProfiler)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/geneInfo.RData")

# Select gene and phenotype of interest
selGene <- "3975" # LHX1
selGene <- "64211" # LHX5
selGene <- "2556" # GABRA3
selGene <- "5630" # PRPH
selPheno <- "Seizure"

selGene <- "64968"
selGene <- "51227"
selGene <- "2114"
selPheno <- "Hypotonia"

selSamples <- metaData_all$ID[metaData_all[,selPheno] == 1]
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
                     Dataset = colnames(logFC_matrix1))
plotDF$Sig <- ifelse((plotDF$Lower > 0) | (plotDF$Upper < 0), "Yes", "No")

# Calculate significance
t <- table(factor(plotDF$Pheno, levels = rev(levels(plotDF$Pheno))), 
           factor(plotDF$Sig, levels = c("No", "Yes")))
fisher.test(t)

# Combine with meta data
plotDF <- inner_join(plotDF, metaData_all, by = c("Dataset" = "ID"))

# Make plot
p <- ggplot(plotDF) +
  geom_point(aes(x = Dataset, y = logFC, 
                 color = logFC)) +
  geom_segment(aes(x = Dataset, xend = Dataset, y = Lower, yend = Upper, 
                   color = logFC)) +
  facet_grid(cols = vars(Pheno), scale = "free", space = "free") +
  scale_color_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                        trans = "pseudo_log") +
  ylab(expression(log[2]~"FC")) +
  xlab("Datasets") +
  ggtitle(name) +
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
# Save plot
ggsave(p, file = "7. PhenoMeta/Gene_Pheno/LHX5.png", width = 6, height = 4)


# Export figure source data
sourceData <- plotDF[,c("Dataset", "logFC", "Lower", "Upper", "Pheno")]
colnames(sourceData) <- c("Dataset", "logFC", "Lower", "Upper", "Phenotype")
write.csv(sourceData, file = "7. PhenoMeta/SourceData_Figure4D.csv",
          row.names = FALSE, quote = FALSE)


# Make an alternative plot
plotDF_alt <- plotDF
plotDF_alt <- arrange(plotDF_alt, by = logFC)
plotDF_alt$Dataset <- factor(plotDF_alt$Dataset, 
                             levels = unique(plotDF_alt$Dataset))

# Set colors
colors <- setNames(c("#BDBDBD", RColorBrewer::brewer.pal(n = 8, name = "Set2")),
                   c(paste0("No ", selPheno),"Intellectual disability", "Hypotonia", 
                     "Global developmental delay",
                     "Microcephaly", "Gait ataxia", "Autism/Autistic Behavior",
                     "Seizure", "Scoliosis"))

# Make plot
p <- ggplot(plotDF_alt) +
  geom_point(aes(x = Dataset, y = logFC, 
                 color = Pheno)) +
  geom_segment(aes(x = Dataset, xend = Dataset, y = Lower, yend = Upper, 
                   color = Pheno)) +
  ylab(expression(log[2]~"FC")) +
  xlab("Datasets") +
  scale_color_manual(values = colors) +
  ggtitle(name) +
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


ggsave(p, file = "7. PhenoMeta/Gene_Pheno/LHX1_alt.png", width = 6, height = 4)

