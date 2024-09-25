# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/GEOData/Analysis")

# Load packages
library(tidyverse)
library(igraph)
library(RCy3)

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

# For each GO term, get the gene that is significant the most times
gene <- rep(NA, length(terms_ordered))
value <- rep(NA, length(terms_ordered))
for (i in 1:length(terms_ordered)){
  selGenes <- GOgenes[[terms_ordered[i]]]
  p_fil <- pvalue_matrix[rownames(pvalue_matrix) %in% selGenes,]
  p_fil[is.na(p_fil)] <- 1
  logFC_fil <- logFC_matrix[rownames(logFC_matrix) %in% selGenes,]
  logFC_fil[is.na(logFC_fil)] <- 0
  SE_fil <- SE_matrix[rownames(SE_matrix) %in% selGenes,]
  SE_fil[is.na(SE_fil)] <- 10
  
  upper <- logFC_fil + 1.96*SE_fil
  lower <- logFC_fil - 1.96*SE_fil
  
  value[i] <- sort(rowSums((lower > 0) | (upper < 0)))[nrow(lower)]
  gene[i] <- names(sort(rowSums((lower > 0) | (upper < 0)))[nrow(lower)])
}

testDF <- data.frame(gene = gene,
                     value = value,
                     GO = terms_ordered)
table(gene)

# Select gene of interest
selGene <- "351"
selGene <- "6622"
selGene <- "3106"
selGene <- "2114"

# Prepare logFC and SE matrix
logFC_matrix1 <- logFC_matrix
SE_matrix1 <- SE_matrix
logFC_matrix1[is.na(logFC_matrix1)] <- 0
SE_matrix1[is.na(SE_matrix1)] <- 0

# Prepare data for plotting
geneName <- geneInfo$Symbol[geneInfo$GeneID == selGene]
plotDF <- data.frame(logFC = logFC_matrix1[selGene,],
                     Upper = logFC_matrix1[selGene,] + 1.96*SE_matrix1[selGene,],
                     Lower = logFC_matrix1[selGene,] - 1.96*SE_matrix1[selGene,],
                     Pvalue = pvalue_matrix[selGene,],
                     StudyID = colnames(logFC_matrix1))
plotDF$Sig <- ifelse((plotDF$Lower > 0) | (plotDF$Upper < 0) , "Yes", "No")
plotDF <- inner_join(plotDF, metaData_all, by = c("StudyID" = "ID"))

# Order the datasets
plotDF$StudyID <- factor(plotDF$StudyID,
                         levels = sets_ordered)

# Make plot
colors <- c("#BDBDBD", "#6A51A3")
p <- ggplot(plotDF) +
  geom_hline(yintercept = 0) +
  geom_segment(aes(x = StudyID, xend = StudyID, y = Lower, yend = Upper, color = Sig)) +
  geom_point(aes(x = StudyID, y = logFC, color = Sig)) +
  ylab(expression(log[2]~FC)) +
  xlab(NULL) +
  ggtitle(geneName) +
  labs(color = NULL) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

# Save plot
ggsave(p, file = paste0("5. OverallMeta/Gene/",geneName, "_logFCs.png"), width = 9, height = 2)


################################################################################

# Network

################################################################################

# Link genes to GO term
networkDF <- NULL
netGenes <- c("351", "6622", "3106", "2114")
GOgenes_fil <- GOgenes[terms_ordered]
for (g in 1:length(netGenes)){
  for (t in 1:length(GOgenes_fil)){
    if (netGenes[g] %in% GOgenes_fil[[t]]){
      networkDF <- rbind.data.frame(networkDF, c(netGenes[g], names(GOgenes_fil)[t]))
    }
  }
}
colnames(networkDF) <- c("geneID", "GOID")

# Add GO name
networkDF <- left_join(networkDF, NES[,1:2], by = c("GOID" = "ID"))
networkDF$Description <- firstup(networkDF$Description)

# Add clusters
networkDF <- left_join(networkDF, clusterDF, by = c("GOID" = "ID"))
networkDF$Cluster[networkDF$Cluster == 1] <- "A"
networkDF$Cluster[networkDF$Cluster == 2] <- "B"
networkDF$Cluster[networkDF$Cluster == 3] <- "A"
networkDF$Cluster[networkDF$Cluster == 4] <- "D"
networkDF$Cluster[networkDF$Cluster == 5] <- "C"

# Add gene symbols
geneInfo$GeneID <- as.character(geneInfo$GeneID)
networkDF <- left_join(networkDF, unique(geneInfo[, c("GeneID", "Symbol")]),
                       by = c("geneID" = "GeneID"))

# Make edges
edges <- networkDF[,c("Symbol", "Description")]
colnames(edges) <- c("from", "to")

# Make nodes
nodes <- unique(data.frame(name = c(networkDF$Symbol, 
                             networkDF$Description),
                    group = c(rep("Gene", nrow(networkDF)),
                              networkDF$Cluster)))

# Make graph
g <- graph_from_data_frame(edges, 
                           directed=FALSE, 
                           vertices=nodes)

# Export graph into Cytoscape
createNetworkFromIgraph(
  g,
  title = "GO network",
  collection = "Genes"
)

