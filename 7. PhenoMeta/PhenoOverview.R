# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/Analysis")

# Load packages
library(tidyverse)

# Load data
load("Data/CleanData/metaData_all.RData")

# Make dataframe for plotting
plotDF <- gather(metaData_all[,13:20])
plotDF$ID <- rep(metaData_all$ID, 8)

# Cluster the studies
ClusMat <- as.matrix(metaData_all[,13:20])
rownames(ClusMat) <- metaData_all$ID
colnames(ClusMat) <- colnames(metaData_all[,13:20])
model_study <- hclust(dist(ClusMat,
                           method = "euclidean"), "ward.D2")
study_ordered <- model_study$labels[model_study$order]
plotDF$ID <- factor(plotDF$ID, levels = study_ordered)

# Add phenotype label
plotDF$value[plotDF$value == 1] <- plotDF$key[plotDF$value == 1]

# Order phenotype
model_pheno <- hclust(dist(t(ClusMat),
                           method = "euclidean"), "ward.D2")
pheno_ordered <- model_pheno$labels[model_pheno$order]
emptyRows <- 3
plotDF <- rbind.data.frame(plotDF, data.frame(key = rep(c(" ", "  ", "   "),each = emptyRows*length(levels(plotDF$ID))),
                                              value = rep(NA,emptyRows*length(levels(plotDF$ID))),
                                              ID = rep(levels(plotDF$ID),emptyRows)))
plotDF$key <- factor(plotDF$key, levels = c(" ", "  ","   ",pheno_ordered))
plotDF$value <- as.character(plotDF$value)
plotDF$value <- factor(plotDF$value, levels = c("0", levels(plotDF$key)[4:11]))

# Set colors
colors <- c("#F0F0F0", RColorBrewer::brewer.pal(n = 8, name = "Set2"))

# Make plot
p <- ggplot(plotDF) +
  geom_tile(aes(y = key, x = ID, fill = value),height = 0.9) +
  scale_fill_manual(values = colors, na.value = "white") +
  xlab("Statistical Comparison") +
  labs(fill = NULL) +
  ylab(NULL) +
  coord_polar() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "right")

# Save plot
ggsave(p, file = "7. PhenoMeta/PhenoOverview/PhenotypeOverview.png", width = 7, height = 4)


