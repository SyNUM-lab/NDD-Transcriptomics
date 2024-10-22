
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(hdf5r)
library(Seurat)
library(glmGamPoi)
library(tidyverse)
library(patchwork)

# Read raw data
expression_matrix <- Read10X_h5("Data/snRNAseq/GSE153184_All_samples_combined_aggr_filtered_feature_bc_matrix.h5")
seurat <- CreateSeuratObject(counts = expression_matrix, 
                             min.cells = 3, min.features = 200)

#==============================================================================#
# Quality control
#==============================================================================#

# Get percentage mitochondrial genes
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

# Visualize QC metrics
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Eliminate low-quality nuclei:
# - Total number of unique molecular identifiers (nCounts_RNA) <500, 
# - Unique detected genes (nFeatures_RNA) <200 
# - Mitochondrial content (percent.mt) >5%
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 5)

#==============================================================================#
# Data processing
#==============================================================================#

# Normalization
seurat <- NormalizeData(seurat)

# Find variable features
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(seurat)

# Scaling
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

# PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
DimPlot(seurat, reduction = "pca") + NoLegend()

# Estimate dimensionality using Elbow plot
ElbowPlot(seurat)

# Cluster
seurat <- FindNeighbors(seurat, dims = 1:6)
seurat <- FindClusters(seurat, resolution = 0.5)

# Perform UMAP
seurat <- RunUMAP(seurat, dims = 1:6)
save(seurat, file = "10. snRNAseq/seurat1.RData")

# Continue with saved data
load("10. snRNAseq/seurat1.RData")

# Extract expression
test <- seurat@assays$RNA$data["Lhx1",]

# Make plot of Lhx1 expression
t <- FeaturePlot(seurat, features = "Lhx1")
plotDF <- t[[1]]$data
p <- ggplot(plotDF) +
  geom_point(aes(x = umap_1, y = umap_2, color = Lhx1),
             size = 0.05) +
  scale_color_gradient(low = "#D9D9D9", high = "#CB181D" ,
                       limits= c(0,2), oob = scales::squish) +
  labs(color = NULL) +
  xlab("UMAP1\n ") +
  ylab("UMAP2") +
  scale_alpha_continuous(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(face = "bold", size= 15),
        axis.title.y = element_text(face = "bold", angle = 90, size = 15))

# save plot
ggsave(p, file = "10. scRNAseq/Lhx1_UMAP.png",
       width = 6, height = 6)

# Make plot of cell type markers
selGenes <- c("Kcnd2", "Gabra6", "Rbfox3",
              "Pcp2", "Calb1", "Car8",
              "Slc24a3", "Esrrg", "Tfap2b",
              "Slc1a3", "Slc1a2", "Aldh1l1",
              "Mbp", "Sox10", "Olig1")

t <- FeaturePlot(seurat, features = selGenes)
markerDF <- data.frame(umap_1 = unlist(lapply(t, function(x) x$data[,1])),
                   umap_2 = unlist(lapply(t, function(x) x$data[,2])),
                   Expr = unlist(lapply(t, function(x) x$data[,4])),
                   Gene = rep(selGenes, each = 74160))


# Combine all plots into single figure
for (i in 1:length(selGenes)){
  sel <- selGenes[i]
  p <- ggplot(markerDF[markerDF$Gene == sel,]) +
    geom_point(aes(x = umap_1, y = umap_2, color = Expr),
               size = 0.05) +
    scale_color_gradient(low = "#D9D9D9", high = "#CB181D" ,
                         limits= c(0,5), oob = scales::squish) +
    ggtitle(sel) +
    labs(color = "Normalized\nexpression") +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  
  if (i == 1){
    finalPlot <- p
  } else{
    finalPlot <- finalPlot + p
    
  }
}
finalPlot <- finalPlot + patchwork::plot_layout(ncol = 3, nrow = 5)

# Save plot
ggsave(finalPlot, file = "10. scRNAseq/QC_UMAP.png",
       width = 10, height = 12)




