# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)

# Load data
load("Data/CleanData/geneInfo.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")

################################################################################
# Down Syndrome
################################################################################

# Get logFCs of chromosome 21 genes
logFC_chr21 <- logFC_matrix[rownames(logFC_matrix) %in% geneInfo$GeneID[which(geneInfo$ChrName == "21")],]

# Formay for plotting
plotDF <- gather(as.data.frame(logFC_chr21))
plotDF$Gene <- rep(rownames(logFC_chr21), ncol(logFC_chr21))
plotDF <- inner_join(plotDF, metaData_all, by = c("key" = "ID"))
plotDF$DS <- ifelse(plotDF$`Disease abbreviation` == "DS", "Down Syndrome", "Other")

# Get IQR
plotDF_median <- plotDF %>% 
  group_by(key) %>%
  reframe(MedianValue = median(value, na.rm = TRUE),
            Upper = quantile(value, 0.75, na.rm = TRUE),
            Lower = quantile(value, 0.25, na.rm =TRUE),
            DS = DS)

plotDF_median <- unique(plotDF_median)
plotDF_median <- arrange(plotDF_median, by = MedianValue)
plotDF_median$key <- factor(plotDF_median$key, levels = unique(plotDF_median$key))

# Set colors
colors <- c("#cf9a02", "#737373")

# Make plot
p <- ggplot(plotDF_median) +
  geom_segment(aes(x = key, xend = key, y = Lower, yend = Upper, color = DS)) +
  geom_point(aes(x = key, y = MedianValue, color = DS), size = 1) +
  xlab("Datasets") +
  ylab(expression(log[2]~"FC")) +
  facet_grid(cols = vars(DS), scale = "free", space = "free") +
  ggtitle("Chromosome 21") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))


# convert to grob
panel_colors <- c("#fdbd07", "#D9D9D9")

gp <- ggplotGrob(p)
for(i in 1:2){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}

# Save plot
ggsave(gp, file = "2. QC/QCmarkers/DS_QC.png", width = 7, height = 5)


################################################################################
# DMD
################################################################################

# Select DMD gene
geneID <- as.character(geneInfo$GeneID[geneInfo$Symbol == "DMD"])
plotDF <- data.frame(logFC = logFC_matrix[geneID,],
                     StudyID = colnames(logFC_matrix))
plotDF <- inner_join(plotDF, metaData_all, by = c("StudyID" = "ID"))
plotDF$DMD <- ifelse(plotDF$`Disease abbreviation` == "DMD", "Duchenne Muscular Dystrophy", "Other")

# Set colors
colors <- c("#7570B3", "#737373")

# Make plot
p <- ggplot(plotDF) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_boxplot(aes(x = DMD, y = logFC, fill = DMD),
               outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(x = DMD, y = logFC, color = DMD), 
             position=position_jitterdodge()) +
  facet_grid(cols = vars(DMD), scale = "free", space = "free") +
  ggtitle("DMD") +
  xlab(NULL) +
  ylab(expression(log[2]~"FC")) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))


# convert to grob
panel_colors <- c("#bfbddc", "#D9D9D9")
gp <- ggplotGrob(p)
for(i in 1:2){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
  
}

# Save plot
ggsave(gp, file = "2. QC/QCmarkers/DMD_QC.png", width = 5, height = 5)


################################################################################
# FXS
################################################################################

# Select FMR1 gene
geneID <- as.character(geneInfo$GeneID[geneInfo$Symbol == "FMR1"])
plotDF <- data.frame(logFC = logFC_matrix[geneID,],
                     StudyID = colnames(logFC_matrix))
plotDF <- inner_join(plotDF, metaData_all, by = c("StudyID" = "ID"))
plotDF$DMD <- ifelse(plotDF$`Disease abbreviation` == "FXS", "Fragile X Syndrome", "Other")

# Set colors
colors <- c("#E7298A", "#737373")

# Make plot
p <- ggplot(plotDF) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_boxplot(aes(x = DMD, y = logFC, fill = DMD),
               outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(x = DMD, y = logFC, color = DMD), 
             position=position_jitterdodge()) +
  facet_grid(cols = vars(DMD), scale = "free", space = "free") +
  ggtitle("FMR1") +
  xlab(NULL) +
  ylab(expression(log[2]~"FC")) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))


# convert to grob
panel_colors <- c("#f07ab6", "#D9D9D9")
gp <- ggplotGrob(p)
for(i in 1:2){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
  
}

# Save plot
ggsave(gp, file = "2. QC/QCmarkers/FXS_QC.png", width = 5, height = 5)


################################################################################
# RTT
################################################################################

# Select MECP2 gene
geneID <- as.character(geneInfo$GeneID[geneInfo$Symbol == "MECP2"])
logFC_fil <- logFC_matrix[geneID,metaData_all$ID[metaData_all$`Disease abbreviation` == "RTT"]]

# Get mutation type for each study
set <- c("GSE107399-1", "GSE107399-2", "GSE107399-3", "GSE113902", "GSE117511",
         "GSE123753-1", "GSE123753-2",  "GSE128380-1", "GSE128380-2", "GSE230714")
mutation <- c("Point mutation", "Point mutation", "Point mutation", "Point mutation",
              "Point mutation", "Deletion", "Deletion","Point mutation", "Point mutation",
              "Point mutation")

# Make dataframe for plotting
plotDF <- data.frame(logFC = logFC_fil,
                     mutation = mutation,
                     set = set)

# Set colors
colors <- c("#D95F02", "#737373")

# Make plot
p <- ggplot(plotDF) +
  geom_bar(aes(x = set, y = logFC, fill = mutation),
           position = position_dodge(), stat = "identity") +
  geom_vline(xintercept = 0, size = 1) +
  facet_grid(cols = vars(mutation), scale = "free", space = "free") +
  xlab("Datsets") +
  ylab(expression(log[2]~FC)) +
  scale_fill_manual(values = colors) +
  ggtitle("MECP2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

# Save plot
ggsave(p, file = "2. QC/QCmarkers/RTT_QC.png", width = 3.3, height = 5)


