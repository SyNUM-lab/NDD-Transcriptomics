# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(patchwork)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Load data
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/metaData_all.RData")
load("5. overallMeta/GSEA/sets_ordered.RData")
load("5. overallMeta/GSEA/terms_ordered.RData")
load("Data/CleanData/geneInfo.RData")

# Get GO annotatation
GOann <- unique(pvalues[,c(1,2)])

# Get all genes from the selected GO terms
genes <- unique(unlist(GOgenes[terms_ordered]))

# Make matrix with GO term-Gene relationship
GOmatrix <- matrix(0, ncol = length(terms_ordered),
                   nrow = length(genes))
rownames(GOmatrix) <- genes
colnames(GOmatrix) <- terms_ordered
for (i in 1:length(terms_ordered)){
  GOmatrix[genes %in% GOgenes[[terms_ordered[i]]],i] <- 1
}

# Cluster the genes
genes_ordered <- rev(genes)

# Make data for plotting
clusterDF$Cluster[clusterDF$Cluster == 1] <- "A"
clusterDF$Cluster[clusterDF$Cluster == 2] <- "B"
clusterDF$Cluster[clusterDF$Cluster == 3] <- "A"
clusterDF$Cluster[clusterDF$Cluster == 4] <- "D"
clusterDF$Cluster[clusterDF$Cluster == 5] <- "C"
GOdf <- NULL
for (i in 1:length(terms_ordered)){
  temp <- data.frame(GeneID = GOgenes[[terms_ordered[i]]],
                     GOID = terms_ordered[i],
                     Cluster = clusterDF$Cluster[clusterDF$ID == terms_ordered[i]])
  
  GOdf <- rbind.data.frame(GOdf, temp)
}
GOdf <- inner_join(GOdf, GOann, by = c("GOID" = "ID"))
GOdf$GeneID <- factor(GOdf$GeneID, levels = genes_ordered)
GOdf$GOID <- factor(GOdf$GOID, levels = terms_ordered)
GOdf$Name <- firstup(GOdf$Description)
GOdf$Name[nchar(GOdf$Name)>50] <- paste0(substring(GOdf$Name[nchar(GOdf$Name)>50],1,47),"...")
GOdf <- arrange(GOdf, by = GOID)
GOdf$Name <- factor(GOdf$Name, levels = unique(GOdf$Name))

# Make plot
colors <- c("#1B9E77","#D95F02", "#E7298A","#E6AB02")
p_GO <- ggplot(GOdf) +
  geom_tile(aes(x = Name, y = GeneID, fill = Cluster)) +
  facet_grid(cols = vars(Cluster), scale = "free", space = "free") +
  theme_bw() +
  scale_fill_manual(values = colors) +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        legend.position = "none")


# Make histogram
p_GOhist <- ggplot(GOdf) +
  geom_bar(aes(x = Name, fill = Cluster)) +
  facet_grid(cols = vars(Cluster), scale = "free", space = "free") +
  ylab("# Genes") +
  theme_bw() +
  scale_fill_manual(values = colors) +
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")


# Prepare gene-logFC data for plotting
logFC_fil <- logFC_matrix[rownames(logFC_matrix) %in% genes,]
logFC_fil[is.na(logFC_fil)] <- 0
plotDF <- gather(as.data.frame(logFC_fil))
plotDF$GeneID <- rep(rownames(logFC_fil), ncol(logFC_fil))
plotDF$key <- factor(plotDF$key, levels = sets_ordered)
plotDF$GeneID <- factor(plotDF$GeneID, levels = genes_ordered)

# Make plot
p_sets <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = GeneID, fill = value)) +
  scale_fill_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                       trans = "pseudo_log", limits = c(-2,2),
                       oob = scales::squish) +
  xlab("Datasets") +
  ylab("Genes") +
  labs(fill = expression(log[2]~"FC")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Combine plots
p_final <- patchwork::plot_spacer() + p_GOhist +
  p_sets + p_GO + patchwork::plot_layout(ncol = 2, nrow = 2,
                                       widths = c(1,1), heights = c(0.2,1))

# Save plot
ggsave(p_final, file = "5. OverallMeta/GSEA/GO_gene_plot.png",
       width = 8, height = 10)

# Write CSV with GO-gene relationship
geneInfo$GeneID <- as.character(geneInfo$GeneID)
GOdf <- left_join(GOdf, unique(geneInfo[,c(1,2)]), by = c("GeneID" = "GeneID"))
printGO <- GOdf[,c("GeneID", "Symbol", "GOID", "Name")]
colnames(printGO) <- c("EntrezGene", "GeneSymbol", "GOID", "GOName")

# Add logFCs to table
logFCdf <- as.data.frame(logFC_matrix)
logFCdf$Gene <- as.character(rownames(logFCdf))
printGO <- inner_join(printGO,logFCdf, by = c("EntrezGene" = "Gene"))

# Write csv
write.csv(printGO, file = "5. OverallMeta/GSEA/GO-gene relationship.csv",
          quote = FALSE, row.names = FALSE, col.names = TRUE)

