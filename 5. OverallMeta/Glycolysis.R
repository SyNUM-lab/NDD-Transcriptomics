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

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Load data
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/geneInfo.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)

# Load WikiPathways
gmt <- clusterProfiler::read.gmt.wp("5. OverallMeta/Glycolysis/wikipathways-20240710-gmt-Homo_sapiens.gmt")
path2gene <- unique(gmt[,c("wpid", "gene")])
path2name <- unique(gmt[,c("wpid", "name")])

# Get the number of significant datasets for each gene
logFC_matrix1 <- logFC_matrix
SE_matrix1 <- SE_matrix
logFC_matrix1[is.na(logFC_matrix1)] <- 0
SE_matrix1[is.na(SE_matrix1)] <- 0
upper = logFC_matrix1 + 1.96*SE_matrix1
lower = logFC_matrix1 - 1.96*SE_matrix1
resultDF <- data.frame(gene = as.character(rownames(pvalue_matrix)),
                       sigValue = rowSums((lower > 0) | (upper < 0)))
resultDF <- inner_join(resultDF, unique(geneInfo[,c(1,2,6)]), by = c("gene" = "GeneID"))
resultDF <- arrange(resultDF, by = desc(sigValue))

# Perform WikiPathways Overrepresentation analysis on top 500 genes
set.seed(123)
WPtest <- enricher(gene = unique(resultDF$gene[1:500]),
                   universe = unique(resultDF$gene),
                   pAdjustMethod = "fdr",
                   pvalueCutoff = Inf,
                   qvalueCutoff = Inf,
                   TERM2GENE = path2gene,
                   TERM2NAME = path2name)

WPresults <- WPtest@result

set.seed(123)
GOtest <- enrichGO(gene = unique(resultDF$gene[1:500]),
                   universe = unique(resultDF$gene),
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = Inf,
                   qvalueCutoff = Inf)

GOresults <- GOtest@result

################################################################################

# Glycolysis

################################################################################

# Select glycolysis genes
selGenes <- path2gene$gene[path2gene$wpid == "WP4628"]
gene_order <- c("SLC2A1", "HK1", "GPI", "PFKM", "ALDOA", "GAPDH",
                "PGK1", "PGAM2", "ENO1", "PKM", "LDHA")

# Prepare dataframe for plotting
plotDF <- resultDF[resultDF$Symbol %in% gene_order,]
plotDF$Symbol <- factor(plotDF$Symbol, levels = gene_order)

# Make colors
binwidth = 1
n_bins <- length(ggplot2:::bin_breaks_width(range(resultDF$sigValue), 
                                            width = binwidth)$breaks) - 1L
colors <- setNames(rev(heat.colors(n_bins)),
                   as.character(1:45))

# Make lollipop diagram
p1 <- ggplot(plotDF) +
  geom_segment(aes(x = Symbol, xend = Symbol, y = 0, yend = sigValue),
               color = "#737373", size = 1) +
  geom_point(aes(x = Symbol, y = sigValue, color = as.character(sigValue)),
             size = 5) +
  coord_cartesian(ylim= c(0,45)) +
  scale_color_manual(values = colors) +
  xlab(NULL) +
  ylab("# significant datasets") +
  theme_bw() +
  theme(legend.position = "none") 

# Save plot
ggsave(p2, file = "5. OverallMeta/Glycolysis/glycolysis1.png", width = 3, height = 5)


# Make histogram
p2 <- ggplot(resultDF) +
  geom_histogram(aes(y = sigValue), fill = rev(heat.colors(n_bins)),
                 binwidth = binwidth, color = "black") +
  coord_cartesian(ylim= c(0,45)) +
  ylab(NULL) +
  xlab("Count") +
  theme_bw()

# Save plot
ggsave(p2, file = "5. OverallMeta/Glycolysis/glycolysis2.png", width = 3, height = 5)
