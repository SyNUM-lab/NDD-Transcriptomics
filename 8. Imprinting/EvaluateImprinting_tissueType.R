# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")
load("Data/CleanData/geneInfo.RData")
load("Data/CleanData/topList.RData")
metaData <- readxl::read_excel("Data/CleanData/MetaData_clean.xlsx")

# Load packages
library(tidyverse)
library(readxl)

# Get imprinted genes: https://www.geneimprint.com/site/genes-by-species
imprintDF <- read_xlsx("Data/ImprintedGenes.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]

# Get imprinted genes
imprinted_genes1 <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
imprinted_genes1 <- imprinted_genes1[!is.na(imprinted_genes1)]
imprinted_genes1 <- unique(as.character(geneInfo$GeneID[geneInfo$Symbol %in% imprinted_genes1]))

up_genes <- list()
down_genes <- list()
both_genes <- list()
all_genes <- list()
for (t in 1:length(topList)){
  testList <- topList[[t]][!is.na(topList[[t]]$`p-value`),]
  up_genes[[t]] <- testList$GeneID[(testList$`p-value` < 0.05) & (testList$log2FC > 0)]
  down_genes[[t]] <- testList$GeneID[(testList$`p-value` < 0.05) & (testList$log2FC < 0)]
  both_genes[[t]] <- testList$GeneID[(testList$`p-value` < 0.05)]
  all_genes[[t]] <- testList$GeneID
}

results_up <- matrix(NA,nrow = length(both_genes), ncol = 4)
results_down <- matrix(NA,nrow = length(both_genes), ncol = 4)
results_both <- matrix(NA,nrow = length(both_genes), ncol = 4)
for (d in 1:length(both_genes)){
  
  # get all genes
  all_genes1 <- all_genes[[d]]
  
  # imprinted genes
  imprinted_genes <- intersect(all_genes1, imprinted_genes1)
  
  # get non-imprinted genes
  nimprinted_genes <- setdiff(all_genes1, imprinted_genes)
  nimprinted_genes <- nimprinted_genes[!is.na(nimprinted_genes)]
  
  # Evaluate upregulated genes:
  
  # get upregulated genes
  up_genes1 <- up_genes[[d]]
  
  # get non-upregulated genes
  nup_genes1 <- setdiff(all_genes1, up_genes1)
  
  
  up_imp <- intersect(up_genes1,imprinted_genes)
  nup_imp <- intersect(nup_genes1,imprinted_genes)
  up_nimp <- intersect(up_genes1, nimprinted_genes)
  nup_nimp <- intersect(nup_genes1, nimprinted_genes)
  
  m <- matrix(c(length(up_imp), 
                length(nup_imp), 
                length(up_nimp),
                length(nup_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_up[d,] <- c(output$p.value,output$estimate,output$conf.int)
  
  # Evaluate downregulated genes:
  
  # get downpregulated genes
  down_genes1 <- down_genes[[d]]
  
  # get non-downregulated genes
  ndown_genes1 <- setdiff(all_genes1, down_genes1)
  
  
  down_imp <- intersect(down_genes1,imprinted_genes)
  ndown_imp <- intersect(ndown_genes1,imprinted_genes)
  down_nimp <- intersect(down_genes1, nimprinted_genes)
  ndown_nimp <- intersect(ndown_genes1, nimprinted_genes)
  
  m <- matrix(c(length(down_imp), 
                length(ndown_imp), 
                length(down_nimp),
                length(ndown_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_down[d,] <- c(output$p.value,output$estimate,output$conf.int)
  
  
  # Evaluate up- and downregulated genes:
  
  # get DEd genes
  both_genes1 <- both_genes[[d]]
  
  # get non-DEd genes
  nboth_genes1 <- setdiff(all_genes1, both_genes1)
  
  
  both_imp <- intersect(both_genes1,imprinted_genes)
  nboth_imp <- intersect(nboth_genes1,imprinted_genes)
  both_nimp <- intersect(both_genes1, nimprinted_genes)
  nboth_nimp <- intersect(nboth_genes1, nimprinted_genes)
  
  m <- matrix(c(length(both_imp), 
                length(nboth_imp), 
                length(both_nimp),
                length(nboth_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_both[d,] <- c(output$p.value,output$estimate,output$conf.int)
}

# Format data
colnames(results_both) <- c("Pvalue", "OR", "Lower", "Upper")
rownames(results_both) <- names(topList)
plotDF_both <- as.data.frame(results_both)
plotDF_both$StudyID <- rownames(results_both)

colnames(results_down) <- c("Pvalue", "OR", "Lower", "Upper")
rownames(results_down) <- names(topList)
plotDF_down <- as.data.frame(results_down)
plotDF_down$StudyID <- rownames(results_down)

colnames(results_up) <- c("Pvalue", "OR", "Lower", "Upper")
rownames(results_up) <- names(topList)
plotDF_up <- as.data.frame(results_up)
plotDF_up$StudyID <- rownames(results_up)

#==============================================================================#
# Make plot
#==============================================================================#

# Odds of differential expression:

# Format data for plotting
plotDF_both <- inner_join(plotDF_both, metaData, by = c("StudyID" = "ID"))
plotDF_both <- plotDF_both[plotDF_both$System != "Primary",]
plotDF_both <- arrange(plotDF_both, by = OR)
plotDF_both$x <- 1:nrow(plotDF_both)

# Load random permutations
load("8. Imprinting/PermResults_invitro.RData")
# Alternative: load("8. Imprinting/PermResults_noninvitro.RData")
plotDF_perm <- gather(as.data.frame(results_both))
plotDF_perm$x <- rep(1:nrow(results_both), ncol(results_both))


# Make plot
p <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_line(data = plotDF_perm, aes(x = x, y = log2(value), group = key), 
            alpha = 0.05, color = "#D94801") +
  geom_point(data = plotDF_both, aes(x = x, y = log2(OR)), 
             alpha = 1, color = "#525252", size = 0.8) +
  geom_line(data = plotDF_both, aes(x = x, y = log2(OR), group = 1),
             alpha = 1, color = "#525252") +
  geom_segment(data = plotDF_both, aes(x = x, xend = x, 
                                       y = log2(Lower), yend = log2(Upper)), 
               alpha = 1, color = "#525252") +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  ylab(expression(log[2]~"OR")) +
  xlab("Datasets") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save plot
ggsave(p, file = "8. Imprinting/Figures/OR_imprinting_invitro.png", width = 7, height = 3)

# Calculate permutation P value
sum(plotDF_both$OR > 1)
sum(colSums(results_both > 1)>=99)