# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

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
load("Data/CleanData/topList.RData")
load("Data/CleanData/geneInfo.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)
geneInfo$ChrName[geneInfo$ChrAcc == "NC_012920.1"] <- "MT"

# Load WikiPathways
geneInfo_fil <- geneInfo[!is.na(geneInfo$ChrName),]
path2gene <- unique(geneInfo[,c("ChrName","GeneID")])

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Perform GSEA for each study
for (i in 1:length(topList)){
  
  # Maked sorted list based on logFC
  gsea_input <- sort(setNames(topList[[i]]$log2FC,topList[[i]]$GeneID), 
                     decreasing = TRUE)
  
  # Perform GSEA
  set.seed(123)
  CHRtest <- GSEA(
    geneList = gsea_input,
    pAdjustMethod = "fdr",
    TERM2GENE = path2gene,
    pvalueCutoff = Inf,
    minGSSize = -Inf,
    maxGSSize = Inf)
  
  # Get results
  CHRresults <- CHRtest@result
  
  # Collect P values and NES
  if (i == 1){
    pvalues <- CHRresults[,c("ID", "pvalue")]
    colnames(pvalues) <- c("ID", names(topList)[i])
    NES <- CHRresults[,c("ID", "NES")]
    colnames(NES) <- c("ID", names(topList)[i])
  } else{
    temp_p <- CHRresults[,c("ID", "pvalue")]
    colnames(temp_p) <- c("ID", names(topList)[i])
    pvalues <- full_join(pvalues, temp_p, by = c("ID" = "ID"))
    
    temp_n <- CHRresults[,c("ID", "NES")]
    colnames(temp_n) <- c("ID", names(topList)[i])
    NES <- full_join(NES, temp_n, by = c("ID" = "ID"))
  }
  
}

save(pvalues, NES, file = "4. GSEA/GSEA_CHR/GSEAresults_CHR.RData")

# Calculate signed P value
pvalues1 <- -log10(pvalues[,-1]) * sign(NES[,-1])
rownames(pvalues1) <- pvalues$ID

# Collect statistics
statDF <- data.frame(
  Chr = pvalues$ID,
  up = rowSums(pvalues1 > -log10(0.05), na.rm = TRUE),
  down = rowSums(pvalues1 < log10(0.05), na.rm = TRUE),
  total = rowSums(abs(pvalues1) > -log10(0.05), na.rm = TRUE))

statDF <- arrange(statDF, by = total)
statDF$Chr <- factor(statDF$Chr, levels = unique(statDF$Chr))


# Make bar plot of # datasets that have P < 0.05:

# Prepare data
barPlot <- data.frame(
  Chr = rep(statDF$Chr,2),
  value = c(statDF$up, statDF$down),
  Sig = c(rep("Upregulated", nrow(statDF)), 
          rep("Downregulated", nrow(statDF)))
)

# Set colors
colors <- setNames(c("#6BAED6", "#BDBDBD", "#FB6A4A"),
                   c("Downregulated", "Unchanged", "Upregulated"))

# Make plot
p1 <- ggplot(barPlot) +
  geom_bar(aes(x = value, y = Chr, fill = Sig),
           stat  ="identity", color ="black") +
  scale_fill_manual(values = colors) +
  ylab("Chromosome") +
  xlab("# P value < 0.05") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(p1, file = "4. GSEA/GSEA_CHR/Chr_bar.png",
       width = 4, height = 4)

plotDF <- gather(pvalues1)
plotDF$Chr <- rep(rownames(pvalues1), ncol(pvalues1))
plotDF$Sig <- "Unchanged"
plotDF$Sig[(plotDF$value < log10(0.05))] <- "Downregulated"
plotDF$Sig[(plotDF$value > -log10(0.05))] <- "Upregulated"
plotDF$Chr <- factor(plotDF$Chr, levels = levels(statDF$Chr))

# colors <- setNames(c("#000072", "#BDBDBD", "#CB181D"),
#                    c("Downregulated", "Unchanged", "Upregulated"))

colors <- setNames(c("#6BAED6", "#BDBDBD", "#FB6A4A"),
                   c("Downregulated", "Unchanged", "Upregulated"))

alphas <- setNames(c(0.6,0.1,0.6),
                   c("Downregulated", "Unchanged", "Upregulated"))
p2 <- ggplot(plotDF) +
  geom_jitter(aes(y = Chr, x = value, color = Sig, alpha = Sig), height = 0.2) +
  scale_color_manual(values = colors) +
  scale_alpha_manual(values = alphas) +
  ylab("Chromosome") +
  xlab(expression("Signed -" ~ log[10]~"P value")) +
  theme_bw() +
  theme(legend.position  ="none")

ggsave(p2, file = "4. GSEA/GSEA_CHR/Chr_point.png",
       width = 6, height = 4)

