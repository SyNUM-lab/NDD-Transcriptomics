# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/Analysis")

# Load packages
library(tidyverse)
library(clusterProfiler)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/geneInfo.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Select phenotype of interest
selPheno <- "Seizure"
selSamples <- metaData_all$ID[metaData_all[,selPheno] == 1]
selGenes <- rownames(pvalue_matrix)[rowSums(pvalue_matrix[,selSamples]<0.05, na.rm = TRUE) > 10]

# Update matrix
logFC_matrix1 <- logFC_matrix
logFC_matrix1[is.na(logFC_matrix1)] <- 0
SE_matrix1 <- SE_matrix
SE_matrix1[is.na(SE_matrix1)] <- 0

# Calculate OR of upregulation, downregulation, and both
stat_both <- matrix(nrow = length(selGenes), ncol = 4)
stat_pos <- matrix(nrow = length(selGenes), ncol = 4)
stat_neg <- matrix(nrow = length(selGenes), ncol = 4)

for (i in 1:length(selGenes)){
  sel <- selGenes[i]
  name <- geneInfo$Symbol[geneInfo$GeneID == sel]
  plotDF <- data.frame(logFC = logFC_matrix1[which(rownames(logFC_matrix1) == sel),],
                       Lower = logFC_matrix1[which(rownames(logFC_matrix1) == sel),] - 1.96*SE_matrix1[which(rownames(logFC_matrix1) == sel),],
                       Upper = logFC_matrix1[which(rownames(logFC_matrix1) == sel),] + 1.96*SE_matrix1[which(rownames(logFC_matrix1) == sel),],
                       Pheno = factor(ifelse(colnames(logFC_matrix1) %in% selSamples, selPheno, paste0("No ", selPheno)),
                                      levels = c(selPheno, paste0("No ", selPheno))),
                       Dataset = colnames(logFC_matrix1))
  
  # Determine significance
  plotDF$Sig_both <- ifelse((plotDF$Lower > 0) | (plotDF$Upper < 0), "Yes", "No")
  plotDF$Sig_pos <- ifelse((plotDF$Lower > 0), "Yes", "No")
  plotDF$Sig_neg <- ifelse((plotDF$Upper < 0), "Yes", "No")
  
  # Fisher's exact test (both up- and downregulation)
  t_both <- fisher.test(table(factor(plotDF$Pheno, levels = rev(levels(plotDF$Pheno))), 
             factor(plotDF$Sig_both, levels = c("No", "Yes"))))
  stat_both[i,] <- c(t_both$estimate, t_both$conf.int[1], 
                     t_both$conf.int[2], t_both$p.value)
  
  # Fisher's exact test (upregulation only)
  t_pos <- fisher.test(table(factor(plotDF$Pheno, levels = rev(levels(plotDF$Pheno))), 
                              factor(plotDF$Sig_pos, levels = c("No", "Yes"))))
  stat_pos[i,] <- c(t_pos$estimate, t_pos$conf.int[1], 
                     t_pos$conf.int[2], t_pos$p.value)
  
  # Fisher's exact test (downregulation only)
  t_neg <- fisher.test(table(factor(plotDF$Pheno, levels = rev(levels(plotDF$Pheno))), 
                              factor(plotDF$Sig_neg, levels = c("No", "Yes"))))
  stat_neg[i,] <- c(t_neg$estimate, t_neg$conf.int[1], 
                     t_neg$conf.int[2], t_neg$p.value)
  
}

# Set column and row names
colNames <- c("Estimate", "Lower", "Upper", "Pvalue")
colnames(stat_both) <- colNames
rownames(stat_both) <- selGenes
colnames(stat_neg) <- colNames
rownames(stat_neg) <- selGenes
colnames(stat_pos) <- colNames
rownames(stat_pos) <- selGenes

# Make data frame
stat_both <- as.data.frame(stat_both)
stat_pos <- as.data.frame(stat_pos)
stat_neg <- as.data.frame(stat_neg)

# Add adj. P value
stat_both$adjP <- p.adjust(stat_both$Pvalue, method = "fdr")
stat_pos$adjP <- p.adjust(stat_pos$Pvalue, method = "fdr")
stat_neg$adjP <- p.adjust(stat_neg$Pvalue, method = "fdr")

# Combine all statistics
stat_all <- data.frame(
  GeneID = selGenes,
  Both = stat_both$Pvalue,
  Neg = stat_neg$Pvalue,
  Pos = stat_pos$Pvalue,
  Both_OR = stat_both$Estimate,
  Neg_OR = stat_neg$Estimate,
  Pos_OR = stat_pos$Estimate
)

stat_all <- inner_join(stat_all, geneInfo[,1:2],
                       by = c("GeneID" = "GeneID"))


# Which genes to include in the plot?
sigGenes <- c(rownames(head(arrange(stat_both[stat_both$Estimate > 1,], by = Pvalue),5)),
              rownames(head(arrange(stat_pos[stat_pos$Estimate > 1,], by = Pvalue),5)),
              rownames(head(arrange(stat_neg[stat_neg$Estimate > 1,], by = Pvalue),6))
              )

# Prepare data for plotting
plotDF <- data.frame(
  GeneID = rep(stat_all$GeneID,3),
  Symbol = rep(stat_all$Symbol, 3),
  Pvalue = c(stat_all$Both, stat_all$Neg, stat_all$Pos),
  OR = c(stat_all$Both_OR, stat_all$Neg_OR, stat_all$Pos_OR),
  Stat = rep(c("Both", "Downregulated", "Upregulated"), each = length(selGenes))
)

plotDF$Stat <- factor(plotDF$Stat, levels = rev(c("Both", "Downregulated", "Upregulated")))
plotDF <- plotDF[plotDF$GeneID %in% sigGenes,]
plotDF$GeneID <- factor(plotDF$GeneID, levels = unique(sigGenes))
plotDF <- arrange(plotDF, by = GeneID)
plotDF$Symbol <- factor(plotDF$Symbol,
                        levels = rev(unique(plotDF$Symbol)))
plotDF$OR[plotDF$OR == Inf] <- 20


# Make plot
p <- ggplot(plotDF) +
  geom_vline(xintercept = -log10(0.05), size = 1, linetype = "dashed",
             color = "#525252") +
  geom_bar(aes(y = Symbol,x = -log10(Pvalue), fill = log2(OR)),
           stat = "Identity", position = position_dodge(),
           color = "black") +
  xlab(expression("-"~log[10]~"P value")) +
  ylab(NULL) +
  labs(fill = expression(log[2]~"OR")) +
  facet_grid(cols = vars(Stat), scale = "free") +
  scale_fill_gradient2(low = "#2171B5", 
                       mid = "white", 
                       high = "#CB181D", 
                       midpoint = 0,
                       limits = c(-3,3),
                       oob = scales::squish) +
  theme_bw()

# Save plot
ggsave(p, file = "7. PhenoMeta/Gene_Pheno/SeizureGenes.png", width = 6, height = 4)


