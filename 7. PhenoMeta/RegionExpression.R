# Data available from:
# https://www.proteinatlas.org/humanproteome/brain/data

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Load packages
library(tidyverse)
library(Seurat)
library(sparseMatrixStats)
library(patchwork)
library(zellkonverter)
library(DelayedArray)
library(SummarizedExperiment)

################################################################################

# scRNA-seq from first trimester developing human  brain
# https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077

################################################################################

# Set data directory
dataDir <- "Data/EarlyDevelopment/"

# Load data
sce <- readH5AD(paste0(dataDir, "f5960e39-50cd-4b31-af54-e70e1b5517a8.h5ad"), 
                use_hdf5 = TRUE)

# Extract LHX1 expression
lhx1 <- assay(sce, "X")["ENSG00000273706", ]

# Extract LHX5 expression
lhx5 <- assay(sce, "X")["ENSG00000089116", ]

# Get annotations/metadata
metaData <- as.data.frame(colData(sce))
all(names(lhx1) == rownames(metaData))

# Add LHX1/5 expression
metaData$LHX1 <- lhx1
metaData$LHX5 <- lhx5

# Save metadata
save(metaData, file = paste0(dataDir, "metaData.RData"))

# Load metadata
load(paste0(dataDir, "metaData.RData"))

#==============================================================================#
# All cells combined
#==============================================================================#

# make copy of metadata
metaData_fil <- metaData

# Prepare metadata for plotting
plotDF <- unique(metaData_fil[,c("Region", "Age")])
plotDF$LHX1_mean <- NA
plotDF$LHX1_n <- NA
plotDF$LHX5_mean <- NA
plotDF$LHX5_n <- NA
for (i in 1:nrow(plotDF)){
  LHX1expr <- (metaData_fil$LHX1 != 0) & (metaData_fil$Region == plotDF$Region[i]) & (metaData_fil$Age == plotDF$Age[i])
  plotDF$LHX1_mean[i] <- mean(metaData_fil$LHX1[which(LHX1expr)])
  plotDF$LHX1_n[i] <- sum(LHX1expr)/length(LHX1expr)
  
  
  LHX5expr <- (metaData_fil$LHX5 != 0) & (metaData_fil$Region == plotDF$Region[i]) & (metaData_fil$Age == plotDF$Age[i])
  plotDF$LHX5_mean[i] <- mean(metaData_fil$LHX5[which(LHX5expr)])
  plotDF$LHX5_n[i] <- sum(LHX5expr)/length(LHX5expr)
}

plotDF$Age_num <- as.numeric(as.character(plotDF$Age))


# Exclude head and brain samples 
plotDF <- plotDF[!(plotDF$Region %in% c("Head", "Brain")),]

# Set order of brain regions
plotDF$Region <- factor(plotDF$Region,
                        levels = c("Hindbrain", "Medulla", "Cerebellum", "Pons",
                                   "Midbrain",
                                   "Forebrain","Diencephalon", "Telencephalon"))

# Set colors
color <- setNames(c(rev(c("#FEE0D2", "#FC9272","#EF3B2C","#A50F15")),
                    "#74C476",
                    rev(c("#C6DBEF", "#6BAED6","#08519C"))),
                  c("Hindbrain", "Medulla", "Cerebellum", "Pons",
                    "Midbrain",
                    "Forebrain","Diencephalon", "Telencephalon"))

# Set shapes
shape <- setNames(c(rep(16,4),
                    15,
                    rep(17,3)),
                  c("Hindbrain", "Medulla", "Cerebellum", "Pons",
                    "Midbrain",
                    "Forebrain","Diencephalon", "Telencephalon"))

# Make plot for LHX1
p <- ggplot(plotDF) +
  geom_point(aes(x = Age_num, y = LHX1_mean, color = Region, size = LHX1_n,
                 shape = Region)) +
  labs(color = "Region",
       size = "Cell\nproportion",
       y = "Mean expression",
       x = "Age (pcw)") +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=1),
         size = guide_legend(order=2)) +
  theme_bw()

# Save plot
ggsave(p, file = "7. PhenoMeta/EarlyDevelopment/LHX1_all.png", width = 7, height = 5)


# Make plot for LHX5
p <- ggplot(plotDF) +
  geom_point(aes(x = Age_num, y = LHX5_mean, color = Region, size = LHX5_n,
                 shape = Region)) +
  labs(color = "Region",
       size = "Cell\nproportion",
       y = "Mean expression",
       x = "Age (pcw)") +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=1),
         size = guide_legend(order=2)) +
  theme_bw()

# Save plot
ggsave(p, file = "7. PhenoMeta/EarlyDevelopment/LHX5_all.png", width = 7, height = 5)


#==============================================================================#
# Neuronal cell fraction only
#==============================================================================#

# filter metadata
metaData_fil <- metaData[metaData$cell_type == "neuron",]

# Prepare metadata for plotting
plotDF <- unique(metaData_fil[,c("Region", "Age")])
plotDF$LHX1_mean <- NA
plotDF$LHX1_n <- NA
plotDF$LHX5_mean <- NA
plotDF$LHX5_n <- NA
for (i in 1:nrow(plotDF)){
  LHX1expr <- (metaData_fil$LHX1 != 0) & (metaData_fil$Region == plotDF$Region[i]) & (metaData_fil$Age == plotDF$Age[i])
  plotDF$LHX1_mean[i] <- mean(metaData_fil$LHX1[which(LHX1expr)])
  plotDF$LHX1_n[i] <- sum(LHX1expr)/length(LHX1expr)
  
  
  LHX5expr <- (metaData_fil$LHX5 != 0) & (metaData_fil$Region == plotDF$Region[i]) & (metaData_fil$Age == plotDF$Age[i])
  plotDF$LHX5_mean[i] <- mean(metaData_fil$LHX5[which(LHX5expr)])
  plotDF$LHX5_n[i] <- sum(LHX5expr)/length(LHX5expr)
}

plotDF$Age_num <- as.numeric(as.character(plotDF$Age))


# Exclude head and brain samples 
plotDF <- plotDF[!(plotDF$Region %in% c("Head", "Brain")),]

# Set order of brain regions
plotDF$Region <- factor(plotDF$Region,
                        levels = c("Hindbrain", "Medulla", "Cerebellum", "Pons",
                                   "Midbrain",
                                   "Forebrain","Diencephalon", "Telencephalon"))

# Set colors
color <- setNames(c(rev(c("#FEE0D2", "#FC9272","#EF3B2C","#A50F15")),
                    "#74C476",
                    rev(c("#C6DBEF", "#6BAED6","#08519C"))),
                  c("Hindbrain", "Medulla", "Cerebellum", "Pons",
                    "Midbrain",
                    "Forebrain","Diencephalon", "Telencephalon"))

# Set shapes
shape <- setNames(c(rep(16,4),
                    15,
                    rep(17,3)),
                  c("Hindbrain", "Medulla", "Cerebellum", "Pons",
                    "Midbrain",
                    "Forebrain","Diencephalon", "Telencephalon"))

# Make plot for LHX1
p <- ggplot(plotDF) +
  geom_point(aes(x = Age_num, y = LHX1_mean, color = Region, size = LHX1_n,
                 shape = Region)) +
  labs(color = "Region",
       size = "Cell\nproportion",
       y = "Mean expression",
       x = "Age (pcw)") +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=1),
         size = guide_legend(order=2)) +
  theme_bw()

# Save plot
ggsave(p, file = "7. PhenoMeta/EarlyDevelopment/LHX1_neuronal.png", width = 7, height = 5)


# Make plot for LHX5
p <- ggplot(plotDF) +
  geom_point(aes(x = Age_num, y = LHX5_mean, color = Region, size = LHX5_n,
                 shape = Region)) +
  labs(color = "Region",
       size = "Cell\nproportion",
       y = "Mean expression",
       x = "Age (pcw)") +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=1),
         size = guide_legend(order=2)) +
  theme_bw()

# Save plot
ggsave(p, file = "7. PhenoMeta/EarlyDevelopment/LHX5_neuronal.png", width = 7, height = 5)


################################################################################

# Brain RNA-seq from Human Protein Atlas

################################################################################

dataDir <- "Data/ProteinAtlas/"

#==============================================================================#
# LHX1
#==============================================================================#

# Select LHX1
selGene <- "ENSG00000273706" # LHX1

dirs <- list.files(dataDir)
outputDF <- NULL
for (d in 1:length(dirs)){
  temp <- data.table::fread(list.files(paste0(dataDir,dirs[d]), full.names = TRUE))
  
  if ("nTPM" %in% colnames(temp)){
    temp <- temp[temp$Gene == selGene, c("Brain region", "nTPM")]
    colnames(temp) <- c("BrainRegion", "Expression")
    temp$Expression <- (temp$Expression - mean(temp$Expression))/sd(temp$Expression)
    temp$Dataset <- dirs[d]
    outputDF <- rbind.data.frame(outputDF, temp)
  }else{
    temp <- temp[temp$Gene == selGene, c("Brain region","Expression energy")]
    colnames(temp) <- c("BrainRegion", "Expression")
    temp$Expression <- (temp$Expression - mean(temp$Expression))/sd(temp$Expression)
    temp$Dataset <- dirs[d]
    outputDF <- rbind.data.frame(outputDF, temp)
  }
}

# Seperate Pons and medulla for plot (only for HPA mouse)
outputDF <- rbind.data.frame(outputDF,
                             data.frame(BrainRegion = "medulla oblongata", 
                               Expression = outputDF[outputDF$BrainRegion == "pons and medulla",c("Expression")],
                               Dataset = outputDF[outputDF$BrainRegion == "pons and medulla",c("Dataset")])
                             )

outputDF[outputDF$BrainRegion == "pons and medulla","BrainRegion"] <- "pons"

# Capitilize first letter for plotting
outputDF$BrainRegion <- firstup(outputDF$BrainRegion)

# Prepare dataset names for plotting
outputDF$Dataset <- sapply(outputDF$Dataset, function(x) {switch(x,
                                                                 "Allen Mouse" = "Allen Mouse\nBrain Atlas",
                                                                 "FANTOM" = "FANTOM",
                                                                 "GTEx" = "GTEx",
                                                                 "HPA Human" = "HPA (Human)",
                                                                 "HPA Mouse" = "HPA (Mouse)",
                                                                 "HPA Pig" = "HPA (Pig)"
)})

outputDF$Dataset <- factor(outputDF$Dataset,
                           levels = rev(c("HPA (Human)",
                                      "FANTOM",
                                      "GTEx",
                                      "HPA (Pig)",
                                      "HPA (Mouse)",
                                      "Allen Mouse\nBrain Atlas")))

# Prepare brain region names for plotting
outputDF$BrainRegion[outputDF$BrainRegion == "Hippocampal formation"] <- "Hippocampal\nformation"
outputDF$BrainRegion[outputDF$BrainRegion == "Medulla oblongata"] <- "Medulla\noblongata"

# Make plot
p <- ggplot(outputDF) +
  geom_tile(aes(y = Dataset, x  = BrainRegion, fill = Expression),
            width = 0.95, height = 0.95) +
  #scale_fill_viridis_c(option = "plasma") +
  scale_fill_gradient(low = "#D9D9D9", high = "#CB181D") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")

ggsave(p, file = "7. PhenoMeta/RegionExpression/LHX1_HPA.png",
     width = 7, height = 3.5)



# Export figure source data
sourceData <- outputDF[,c("Dataset", "BrainRegion", "Expression")]
sourceData$Dataset <- str_remove(sourceData$Dataset, "\n")
sourceData$BrainRegion <- str_remove(sourceData$BrainRegion, "\n")
write.csv(sourceData, file = "7. PhenoMeta/SourceData_Figure4E.csv",
          row.names = FALSE, quote = FALSE)

#==============================================================================#
# LHX5
#==============================================================================#

# Select LHX5
selGene <- "ENSG00000089116" # LHX5

dirs <- list.files(dataDir)
outputDF <- NULL
for (d in 1:length(dirs)){
  temp <- data.table::fread(list.files(paste0(dataDir,dirs[d]), full.names = TRUE))
  
  if ("nTPM" %in% colnames(temp)){
    temp <- temp[temp$Gene == selGene, c("Brain region", "nTPM")]
    colnames(temp) <- c("BrainRegion", "Expression")
    temp$Expression <- (temp$Expression - mean(temp$Expression))/sd(temp$Expression)
    temp$Dataset <- dirs[d]
    outputDF <- rbind.data.frame(outputDF, temp)
  }else{
    temp <- temp[temp$Gene == selGene, c("Brain region","Expression energy")]
    colnames(temp) <- c("BrainRegion", "Expression")
    temp$Expression <- (temp$Expression - mean(temp$Expression))/sd(temp$Expression)
    temp$Dataset <- dirs[d]
    outputDF <- rbind.data.frame(outputDF, temp)
  }
}

# Seperate Pons and medulla for plot (only for HPA mouse)
outputDF <- rbind.data.frame(outputDF,
                             data.frame(BrainRegion = "medulla oblongata", 
                                        Expression = outputDF[outputDF$BrainRegion == "pons and medulla",c("Expression")],
                                        Dataset = outputDF[outputDF$BrainRegion == "pons and medulla",c("Dataset")])
)

outputDF[outputDF$BrainRegion == "pons and medulla","BrainRegion"] <- "pons"

# Capitilize first letter for plotting
outputDF$BrainRegion <- firstup(outputDF$BrainRegion)

# Prepare dataset names for plotting
outputDF$Dataset <- sapply(outputDF$Dataset, function(x) {switch(x,
                                                                 "Allen Mouse" = "Allen Mouse\nBrain Atlas",
                                                                 "FANTOM" = "FANTOM",
                                                                 "GTEx" = "GTEx",
                                                                 "HPA Human" = "HPA (Human)",
                                                                 "HPA Mouse" = "HPA (Mouse)",
                                                                 "HPA Pig" = "HPA (Pig)"
)})

outputDF$Dataset <- factor(outputDF$Dataset,
                           levels = rev(c("HPA (Human)",
                                          "FANTOM",
                                          "GTEx",
                                          "HPA (Pig)",
                                          "HPA (Mouse)",
                                          "Allen Mouse\nBrain Atlas")))

# Prepare brain region names for plotting
outputDF$BrainRegion[outputDF$BrainRegion == "Hippocampal formation"] <- "Hippocampal\nformation"
outputDF$BrainRegion[outputDF$BrainRegion == "Medulla oblongata"] <- "Medulla\noblongata"

# Make plot
p <- ggplot(outputDF) +
  geom_tile(aes(y = Dataset, x  = BrainRegion, fill = Expression),
            width = 0.95, height = 0.95) +
  #scale_fill_viridis_c(option = "plasma") +
  scale_fill_gradient(low = "#D9D9D9", high = "#CB181D") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")

ggsave(p, file = "7. PhenoMeta/RegionExpression/LHX5_HPA.png",
       width = 7, height = 3.5)


################################################################################

# Cerebellar Vermis from Human brain Cell Atlas v1.0

################################################################################

# Load scRNA-seq data
exprData <- read_rds("Data/Cerebellum/Cerebellar Vermis.rds")

#==============================================================================#
# LHX1 expression
#==============================================================================#

# Plot LHX1 expression
LHX1 <- FeaturePlot(exprData, features = "ENSG00000273706")

plotDF <- LHX1[[1]]$data
colnames(plotDF) <- c("UMAP_1", "UMAP_2", "ident", "Gene")

p <- ggplot(plotDF) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = Gene),
             size = 0.05) +
  scale_color_gradient(low = "#D9D9D9", high = "#CB181D" ,
                       limits= c(0,5), oob = scales::squish) +
  labs(color = NULL) +
  xlab("UMAP1\n ") +
  ylab("UMAP2") +
  scale_alpha_continuous(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(face = "bold", size= 15),
        axis.title.y = element_text(face = "bold", angle = 90, size = 15))

ggsave(p, file = "7. PhenoMeta/RegionExpression/LHX1_UMAP_CBV_HBA.png",
       width = 6, height = 6)


# Export figure source data
sourceData <- plotDF[,c("UMAP_1", "UMAP_2", "Gene")]
colnames(sourceData) <- c("UMAP1", "UMAP2", "Expr")
write.csv(sourceData, file = "7. PhenoMeta/SourceData_Figure4F.csv",
          row.names = FALSE, quote = FALSE)

#==============================================================================#
# Purkinje cell markers expression: RORA, HS6ST3
#==============================================================================#

# Load gene information
load("Data/CleanData/geneInfo.RData")

# Get expression
selGene <- FeaturePlot(exprData, features = "ENSG00000069667") # RORA
selGene <- FeaturePlot(exprData, features = "ENSG00000185352") # HS6ST3
selGene <- FeaturePlot(exprData, features = "ENSG00000158321") # AUTS2
selGene <- FeaturePlot(exprData, features = "ENSG00000117114") # LPHN2

# Prepare data for plotting
plotDF <- selGene[[1]]$data
colnames(plotDF) <- c("UMAP_1", "UMAP_2", "ident", "Gene")

# Make plot
p <- ggplot(plotDF) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = Gene),
             size = 0.05) +
  scale_color_gradient(low = "#D9D9D9", high = "#CB181D",
                       limits= c(0,3), oob = scales::squish) +
  labs(color = NULL) +
  xlab("UMAP1\n ") +
  ylab("UMAP2") +
  scale_alpha_continuous(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(face = "bold", size= 15),
        axis.title.y = element_text(face = "bold", angle = 90, size = 15))

# save plot
ggsave(p, file = "7. PhenoMeta/RegionExpression/LPHN2_UMAP_CBV_HBA.png",
       width = 6, height = 6)



# Cell type
plotDF$CellType <- exprData$cell_type
p <- ggplot(plotDF) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = CellType),
             size = 0.05) +
  labs(color = NULL) +
  xlab("UMAP1\n ") +
  ylab("UMAP2") +
  scale_alpha_continuous(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(face = "bold", size= 15),
        axis.title.y = element_text(face = "bold", angle = 90, size = 15))


# Supercluster
plotDF$Supercluster <- exprData$supercluster_term
p <- ggplot(plotDF) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = Supercluster),
             size = 0.05) +
  labs(color = NULL) +
  xlab("UMAP1\n ") +
  ylab("UMAP2") +
  scale_alpha_continuous(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(face = "bold", size= 15),
        axis.title.y = element_text(face = "bold", angle = 90, size = 15))


# Percent of neurons expressed in
sum(plotDF$Gene[plotDF$CellType == "neuron"] > 0)/length(plotDF$Gene[plotDF$CellType == "neuron"])
sum(plotDF$Gene[plotDF$Supercluster == "Cerebellar inhibitory"] > 0)/length(plotDF$Gene[plotDF$Supercluster == "Cerebellar inhibitory"])
sum(plotDF$Gene[plotDF$Supercluster == "Upper rhombic lip"] > 0)/length(plotDF$Gene[plotDF$Supercluster == "Upper rhombic lip"])


