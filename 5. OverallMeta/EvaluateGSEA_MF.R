# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO_MF.RData")

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Prepare meta data:

# 1) Add disease info
selDis <- c("ASD", "FXS", "DMD", "DS", "RTT")
metaData_all$ASD <- ifelse(metaData_all$`Disease abbreviation` == "ASD",1,0)
metaData_all$FXS <- ifelse(metaData_all$`Disease abbreviation` == "FXS",1,0)
metaData_all$DMD <- ifelse(metaData_all$`Disease abbreviation` == "DMD",1,0)
metaData_all$DS <- ifelse(metaData_all$`Disease abbreviation` == "DS",1,0)
metaData_all$RTT <- ifelse(metaData_all$`Disease abbreviation` == "RTT",1,0)

# 2) Add tissue info
neural_tissue <- c("Astrocyte", "BMEC; Astrocyte; Neuron", "Brain organoid",
                   "Cerebral organoid", "Cortical organoid", "Cortical progenitor cell", 
                   "Early neuroectoderm", "iPSC; Neuroepithelium; Neuron", "Late neuroectoderm",
                   "Neuron", "NPC", "Telencephalic organoid", "Brain", "Cerebral granule cell",
                   "Neuroepithelial stem cell"
)
immune_tissue <- c("White blood cell", "PBMC", "Leukocyte", "Lymphocyte", 
                   "Lymphoblastoid cell", "B cell",
                   "Blood", "Microglia cell", "Microglia-like cell")

metaData_all$Neural <- ifelse(metaData_all$Tissue %in% neural_tissue,1,0)
metaData_all$Immune <- ifelse(metaData_all$Tissue %in% immune_tissue,1,0)
metaData_all <- arrange(metaData_all, by = ID)
all(colnames(pvalues[,3:153]) == metaData_all$ID)


# Format p values (GSEA)
p_matrix <- as.matrix(pvalues[,3:153])
p_matrix[is.na(p_matrix)] <- 1
rownames(p_matrix) <- pvalues$ID

# calculate adjusted p-values (GSEA)
adjp_matrix <- matrix(p.adjust(p_matrix, method = "fdr"),
                      nrow = nrow(p_matrix), ncol = ncol(p_matrix))
rownames(adjp_matrix) <- rownames(p_matrix)
colnames(adjp_matrix) <- colnames(p_matrix)

scores <- sort(setNames(rowSums(adjp_matrix < 0.05), rownames(adjp_matrix)))
rownames(pvalues) <- pvalues$ID
rownames(NES) <- NES$ID

cbind(pvalues[rev(names(tail(scores,30))), 1:2],rev(tail(scores,30)))

#==============================================================================#
# Make similarity matrix of GO terms
simMatrix <- calculateSimMatrix(
  rownames(adjp_matrix),
  orgdb="org.Hs.eg.db",
  ont = "MF",
  method = "Resnik"
)

# Select GO term with highest number of significant values
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")

save(reducedTerms, file = "5. overallMeta/GSEA/reducedTerms_MF.RData")
#==============================================================================#

load("5. OverallMeta/GSEA/reducedTerms_MF.RData")
selTerms <- reducedTerms[,c("parent", "parentTerm", "score")]
selTerms <- selTerms %>%
  group_by(parent) %>%
  mutate(score = max(score))
rownames(selTerms) <- NULL
selTerms <- unique(selTerms)

# Plot top 30 terms (with most # of significance)
plotTerms <- head(selTerms$parent,20)
pvalues_plot <- pvalues[plotTerms,]
NES_plot <- NES[plotTerms,]

# Prepare data for plotting
plotDF <- gather(pvalues_plot[,3:153])
plotDF$ID <- rep(pvalues_plot$ID,ncol(pvalues_plot)-2)
plotDF$Name <- rep(pvalues_plot$Description,ncol(pvalues_plot)-2)
plotDF$NES <- gather(NES_plot[,3:153])$value
plotDF$value[is.na(plotDF$value)] <- 1
plotDF$NES[is.na(plotDF$NES)] <- 0

# Cluster the GO terms
values_p <- -log10(pvalues_plot[,3:153]) * sign(NES_plot[,3:153])
values_p[is.na(values_p)] <- 0
values_p[values_p > 5] <- 5
values_p[values_p < -5] <- -5
model <- hclust(dist(values_p), "ward.D2")
terms_ordered <- model$labels[model$order]
clusterDF <- data.frame(ID = names(cutree(model, k = 5)),
                        Cluster = cutree(model, k = 5))
save(terms_ordered, clusterDF, file = "5. overallMeta/GSEA/terms_ordered_MF.RData")

# Put the terms and datasets in correct order
load("5. overallMeta/GSEA/sets_ordered.RData")
plotDF <- inner_join(plotDF, clusterDF, by = c("ID" = "ID"))
plotDF$key <- factor(plotDF$key, levels = sets_ordered)
plotDF$ID <- factor(plotDF$ID, levels = terms_ordered)

plotDF$Name <- firstup(plotDF$Name)
plotDF <- plotDF[!is.na(plotDF$Name),]
plotDF$Name[nchar(plotDF$Name)>50] <- paste0(substring(plotDF$Name[nchar(plotDF$Name)>50],1,47),"...")

plotDF <- arrange(plotDF, by = ID)
plotDF$Name <- factor(plotDF$Name, levels = unique(plotDF$Name))
plotDF$Cluster[plotDF$Cluster == 1] <- "B"
plotDF$Cluster[plotDF$Cluster == 2] <- "A"
plotDF$Cluster[plotDF$Cluster == 3] <- "D"
plotDF$Cluster[plotDF$Cluster == 4] <- "C"
plotDF$Cluster[plotDF$Cluster == 5] <- " "
plotDF$Cluster <- factor(plotDF$Cluster,
                         levels = c("A", "B", "C", "D", " "))

# NOTE: cluster 1 and 3 were combined given their similariy in logFCs and 
# included genes

# Make plot
p <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log10(value)*sign(NES)), height = 0.9) +
  facet_grid(rows = vars(Cluster), scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                       trans = "pseudo_log") +
  ylab(NULL) +
  xlab("Datasets") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

panel_colors <- c("#1B9E77","#D95F02", "#E7298A","#E6AB02", "lightgrey")

# convert to grob
gp <- ggplotGrob(p)
for(i in 1:4){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}

# Save plot
ggsave(gp, file = "5. overallMeta/GSEA/GSEAplot_all_MF.png", width = 10, height = 5)


# Get legends:
legendPlot <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log10(value)*sign(NES)), height = 0.9) +
  facet_grid(rows = vars(Cluster), scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                       trans = "pseudo_log") +
  ylab(NULL) +
  xlab("Datasets") +
  theme_bw() +
  theme(legend.position = "right",
        strip.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "5. overallMeta/GSEA/GSEAplot_all_MF_legend.png", width = 4, height = 4)


