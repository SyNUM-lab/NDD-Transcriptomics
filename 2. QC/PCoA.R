# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/statistics_matrix.RData")
logFC_matrix[is.na(logFC_matrix)] <- 0
pvalue_matrix[is.na(pvalue_matrix)] <- 1


# Calculate pairwise correlations (rank-based, spearman)
corMatrix <- matrix(NA,ncol(pvalue_matrix), ncol(pvalue_matrix))
for (i in 1:ncol(pvalue_matrix)){
  corMatrix[,i] <- apply(pvalue_matrix,2,function(x) {cor(x, pvalue_matrix[,i], 
                                                          method = "spearman")})
}
colnames(corMatrix) <- colnames(pvalue_matrix)
rownames(corMatrix) <- colnames(pvalue_matrix)
save(corMatrix, file = "2. QC/PCoA/corMatrix.RData")


################################################################################

# PCoA

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)

# Load data
load("2. QC/PCoA/corMatrix.RData")
load("Data/CleanData/metaData_all.RData")
geneInfo <- data.table::fread("Data/CleanData/Human.GRCh38.p13.annot.tsv")
load("Data/CleanData/statistics_matrix.RData")
logFC_matrix[is.na(logFC_matrix)] <- 0
pvalue_matrix[is.na(pvalue_matrix)] <- 1

# Convert correlations to distance  
corMatrix1 <- 1-corMatrix

# Scale distance matrix
corMatrix_scaled <- t(t(corMatrix1 - rowMeans(corMatrix1)) - colMeans(corMatrix1)) + mean(corMatrix1)

# Perform PCA
pca <- prcomp(corMatrix_scaled, 
              retx = TRUE, # Give rotated variables (PCA loadings)
              center = FALSE,
              scale = FALSE) # Variables are scaled to unit variance

# Get expained variances
expl_var <- round((pca$sdev^2)/sum(pca$sdev^2),3)

# Combine sample info in PCA scores
plotPCA <- as.data.frame(pca$x)[1:5]
plotPCA$studyID <- rownames(plotPCA)
plotPCA <- inner_join(plotPCA, metaData_all, by = c("studyID" = "ID"))


# Neural cells
neural_tissue <- c("Astrocyte", "BMEC; Astrocyte; Neuron", "Brain organoid",
                   "Cerebral organoid", "Cortical organoid", "Cortical progenitor cell", 
                   "Early neuroectoderm", "iPSC; Neuroepithelium; Neuron", "Late neuroectoderm",
                   "Neuron", "NPC", "Telencephalic organoid", "Brain", "Cerebral granule cell",
                   "Neuroepithelial stem cell"
)

# Immune cells
immune_tissue <- c("White blood cell", "PBMC", "Leukocyte", "Lymphocyte", 
                  "Lymphoblastoid cell", "B cell",
                  "Blood", "Microglia cell", "Microglia-like cell")


# Set colors
plotPCA$color <- ifelse(plotPCA$Tissue %in% immune_tissue, "Immune tissue/cells", "Other")
plotPCA$color[plotPCA$Tissue %in% neural_tissue] <- "Neural tissue/cells"
colors <- setNames(c("#CB181D", "#2171B5", "#969696"),
                   c("Immune tissue/cells", "Neural tissue/cells", "Other"))

# Make plot
p <- ggplot(plotPCA) +
  geom_point(aes(x = PC1, y = PC2, color = color),
             size = 2.5) +
  theme_void() +
  labs(color = NULL) +
  xlab(paste0("PCoA1 (", expl_var[1]*100, "%)")) +
  ylab(paste0("PCoA2 (", expl_var[2]*100, "%)")) +
  scale_color_manual(values = colors) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12,
                                    angle = 90))


# Save plot
ggsave(p, file = "2. QC/PCoA/PCoA.png", width = 6, height = 5)


################################################################################

# Plot individual genes

################################################################################

# Get gene P value ranks
pvalue_rank <- matrix(NA, nrow = nrow(pvalue_matrix), ncol = ncol(pvalue_matrix))
for (c in 1:ncol(pvalue_matrix)){
  pvalue_rank[,c] <- rank(-1*pvalue_matrix[,c])
}
rownames(pvalue_rank) <- rownames(pvalue_matrix)
colnames(pvalue_rank) <- colnames(pvalue_matrix)


# Test for association with immune/neural cells
pvalue_wc <- rep(NA, nrow(pvalue_rank))
statistic_wc <- rep(NA, nrow(pvalue_rank))

for (r in 1:nrow(pvalue_rank)){
  pvalue_wc[r] <- wilcox.test(pvalue_rank[r,plotPCA$studyID[plotPCA$color == "Immune tissue/cells"]],
                             pvalue_rank[r,plotPCA$studyID[plotPCA$color == "Neural tissue/cells"]])$p.value
  
  statistic_wc[r] <- median(pvalue_rank[r,plotPCA$studyID[plotPCA$color == "Immune tissue/cells"]]) - median( pvalue_rank[r,plotPCA$studyID[plotPCA$color == "Neural tissue/cells"]])
}
names(pvalue_wc) <- rownames(pvalue_rank)
names(statistic_wc) <- rownames(pvalue_rank)


# Find gene that is associated with immune cells
head(sort(pvalue_wc[statistic_wc > 0]),20)
selGene <- "6977"
plotPCA$TRGV4 <- pvalue_rank[selGene,]

# Make plot
p <- ggplot(plotPCA) +
  geom_point(aes(x = PC1, y = PC2, color = TRGV4),
             size = 2.5) +
  theme_void() +
  scale_color_gradient(low = "#F0F0F0", high = "#EF3B2C") +
  ggtitle(geneInfo$Symbol[geneInfo$GeneID == selGene][1]) +
  labs(color = "P value\nRank") +
  xlab(paste0("PCoA1 (", expl_var[1]*100, "%)")) +
  ylab(paste0("PCoA2 (", expl_var[2]*100, "%)")) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12,
                                    angle = 90))

# Save plot
ggsave(p, file = "2. QC/PCoA/PCoA_immune.png", width = 6, height = 5)


# Find gene that is associated with neural cells
head(sort(pvalue_wc[statistic_wc < 0]),20)
selGene <- "1137"
plotPCA$CHRNA4 <- pvalue_rank[selGene,]

# Make plot
p <- ggplot(plotPCA) +
  geom_point(aes(x = PC1, y = PC2, color = CHRNA4),
             size = 2.5) +
  theme_void() +
  scale_color_gradient(low = "#F0F0F0", high = "#4292C6") +
  ggtitle(geneInfo$Symbol[geneInfo$GeneID == selGene][1]) +
  labs(color = "P value\nRank") +
  xlab(paste0("PCoA1 (", expl_var[1]*100, "%)")) +
  ylab(paste0("PCoA2 (", expl_var[2]*100, "%)")) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12,
                                    angle = 90))

# Save plot
ggsave(p, file = "2. QC/PCoA/PCoA_neural.png", width = 6, height = 5)

# Export figure source data
sourceData <- plotPCA[,c("studyID","PC1", "PC2", "color", "TRGV4", "CHRNA4")]
colnames(sourceData) <- c("Dataset","PCoA1", "PCoA2", "Tissue", "TRGV4", "CHRNA4")
write.csv(sourceData, file = "2. QC/SourceData_Figure1D.csv",
          row.names = FALSE, quote = FALSE)

################################################################################

# Hierarchical clustering

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(ggdendro)

# Load data
load("2. QC/PCoA/corMatrix.RData")
load("Data/CleanData/metaData_all.RData")
geneInfo <- data.table::fread("Data/CleanData/Human.GRCh38.p13.annot.tsv")
load("Data/CleanData/statistics_matrix.RData")
logFC_matrix[is.na(logFC_matrix)] <- 0
pvalue_matrix[is.na(pvalue_matrix)] <- 1

# First, set cell type as color...

# Neural cells
neural_tissue <- c("Astrocyte", "BMEC; Astrocyte; Neuron", "Brain organoid",
                   "Cerebral organoid", "Cortical organoid", "Cortical progenitor cell", 
                   "Early neuroectoderm", "iPSC; Neuroepithelium; Neuron", "Late neuroectoderm",
                   "Neuron", "NPC", "Telencephalic organoid", "Brain", "Cerebral granule cell",
                   "Neuroepithelial stem cell"
)

# Immune cells
immune_tissue <- c("White blood cell", "PBMC", "Leukocyte", "Lymphocyte", 
                   "Lymphoblastoid cell", "B cell",
                   "Blood", "Microglia cell", "Microglia-like cell")


# Set colors
metaData_all$color <- ifelse(metaData_all$Tissue %in% immune_tissue, "Immune tissue/cells", "Other")
metaData_all$color[metaData_all$Tissue %in% neural_tissue] <- "Neural tissue/cells"

# Second, set the disease as color...

metaData_all$color2 <- metaData_all$Disease
metaData_all$color2[!(metaData_all$color2 %in% c("Down Syndrome",
                                               "Rett Syndrome",
                                               "Fragile X Syndrome",
                                               "Duchenne Muscular Dystrophy"))] <- "Other"

# Perform hierarchical clustering
model <- hclust(as.dist(1-corMatrix), method = "ward.D2")

# Make dendrogram
dend <- as.dendrogram(model)
dend_data <- dendro_data(dend, type = "rectangle")
dend_line <- dend_data$segments
dend_point <- dend_data$labels
dend_point <- inner_join(dend_point,metaData_all, by = c("label" = "ID"))

# Set color values
colors <- setNames(c("#CB181D", "#2171B5", "white",
                     "#D95F02", "#7570B3", "#E7298A", "#E6AB02"),
                   c("Immune tissue/cells", "Neural tissue/cells", "Other", 
                     "Rett Syndrome", "Duchenne Muscular Dystrophy", 
                     "Fragile X Syndrome", "Down Syndrome"))

# Set order of colors
dend_point$color <- factor(dend_point$color,
                           levels = c("Immune tissue/cells", "Neural tissue/cells",
                                      "Duchenne Muscular Dystrophy", "Down Syndrome",
                                      "Fragile X Syndrome", "Rett Syndrome", 
                                      "Other"))
dend_point$color2 <- factor(dend_point$color2,
                           levels = c("Immune tissue/cells", "Neural tissue/cells",
                                      "Duchenne Muscular Dystrophy", "Down Syndrome",
                                      "Fragile X Syndrome", "Rett Syndrome", 
                                      "Other"))

# Make plot
p <- ggplot() + 
  geom_segment(data = dend_line, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_tile(data = dend_point, aes(x, y-0.2, fill = color), 
            height = 0.3, width = 1) +
  geom_tile(data = dend_point, aes(x, y-0.55, fill = color2), 
            height = 0.3, width = 1) +
  scale_fill_manual(values  = colors) +
  labs(fill = NULL) +
  theme_void()

# Save plot
ggsave(p, file = "2. QC/PCoA/dendrogram.png", width = 8, height = 5)
