# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(meta)

# Load data
load("Data/CleanData/geneInfo.RData")
load("Data/CleanData/statistics_matrix.RData")
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")

################################################################################

# Meta analysis: repeat for DS, DMD, RTT, and FX

################################################################################


# Select statistics
disease <- "DS"
logFC_matrix_sel <- logFC_matrix[,metaData_all$ID[metaData_all$`Disease abbreviation` == disease]]
SE_matrix_sel <- SE_matrix[,metaData_all$ID[metaData_all$`Disease abbreviation` == disease]]
pvalue_matrix_sel <- pvalue_matrix[,metaData_all$ID[metaData_all$`Disease abbreviation` == disease]]

# Remove genes with NA statistics
keep <- rownames(logFC_matrix_sel)[(rowSums(is.na(logFC_matrix_sel))==0) &
                                     (rowSums(is.na(SE_matrix_sel))==0) &
                                     (rowSums(is.na(pvalue_matrix_sel))==0)]

logFC_matrix_sel <- logFC_matrix_sel[keep,]
SE_matrix_sel <- SE_matrix_sel[keep,]
pvalue_matrix_sel <- pvalue_matrix_sel[keep,]

# Perform meta analysis
output <- NULL
for (i in 1:nrow(logFC_matrix_sel)){
  m.gen <- tryCatch({
    metagen(TE = logFC_matrix_sel[i,],
            seTE = SE_matrix_sel[i,],
            studlab = colnames(logFC_matrix_sel),
            sm = "MD",
            fixed = FALSE,
            random = TRUE,
            hakn = TRUE)
  }, error = function(cond){
    NULL
  })
  
  if (!is.null(m.gen)){
    output <- rbind(output,
                    c(rownames(logFC_matrix_sel)[i],
                      m.gen$pval.random, 
                      m.gen$statistic.random, 
                      m.gen$lower.random, 
                      m.gen$upper.random))
  }
}
colnames(output) <- c("Gene", "Pvalue", "Statistic", "Lower", "Upper")
outputDF <- as.data.frame(output)
save(outputDF, file = "6. DiseaseMeta/Gene/outputDF_DS.RData")



################################################################################

# Make plot: repeat for DS, DMD, RTT, and FX

################################################################################

disease <- "DMD"

# Load data
load(paste0("6. DiseaseMeta/Gene/outputDF_",disease,".RData"))

# Prepare data
outputDF <- outputDF[!is.na(outputDF$Pvalue),]
outputDF$Gene <- as.character(outputDF$Gene)
geneInfo$GeneID <- as.character(geneInfo$GeneID)
outputDF <- left_join(outputDF, geneInfo, by = c("Gene" = "GeneID"))

# For DS: exclude chromosome 21 genes
#outputDF <- outputDF[outputDF$ChrName != "21",]

# Sort based on P values
outputDF <- arrange(outputDF, by = Pvalue)

# Select gene
selGene <- "5217"

# Select logFCs and SEs
logFC_matrix[is.na(logFC_matrix)] <- 0
SE_matrix[is.na(SE_matrix)] <- 0
plotDF <- data.frame(logFC = logFC_matrix[selGene,],
                     Upper = logFC_matrix[selGene,] + 1.96*SE_matrix[selGene,],
                     Lower = logFC_matrix[selGene,] - 1.96*SE_matrix[selGene,],
                     StudyID = names(logFC_matrix[selGene,]))

# Prepare data for plotting
plotDF <- inner_join(plotDF, metaData_all, by = c("StudyID" = "ID"))
plotDF$`Disease abbreviation`[!(plotDF$`Disease abbreviation` %in% c("RTT", "FXS",
                                                                     "DMD", "DS"))] <- "Other"
plotDF$Disease[plotDF$`Disease abbreviation` == "Other"] <- "Other"
plotDF$Disease1 <- ifelse(plotDF$`Disease abbreviation` == disease, "a","b")
plotDF$StudyID <- factor(plotDF$StudyID)

# Set colors
diseaseColors <- setNames(c("#D95F02", "#7570B3", "#E7298A", "#E6AB02"),
                          c("RTT", "DMD", "FXS", "DS"))
colors <- c(diseaseColors[disease], "#737373")
names(colors) <- c("a", "b")
plotDF$BG <- ifelse(plotDF$`Disease abbreviation`==disease,Inf,0)

# Make plot
p <- ggplot(plotDF[plotDF$Disease != "Other",]) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, 
            xmin = -BG, xmax = BG), 
                   fill = colors[1], alpha = 0.005) +
  geom_vline(xintercept = 0, size = 1, color = "#D9D9D9") +
  geom_segment(aes(y = StudyID, yend = StudyID, x = Lower, xend = Upper, 
                   color = Disease1)) +
  geom_point(aes(y = StudyID, x = logFC, color = Disease1)) +
  facet_grid(rows = vars(`Disease abbreviation`), space = "free", scale = "free") +
  scale_color_manual(values = colors) +
  xlab(expression(log[2]~"FC")) +
  ylab("Datasets") +
  ggtitle(geneInfo$Symbol[geneInfo$GeneID == selGene]) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

# Save plot
ggsave(p, file = paste0("6. DiseaseMeta/Gene/", disease, ".png"), 
       width = 3, height = 7)