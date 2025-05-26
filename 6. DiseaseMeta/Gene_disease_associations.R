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
disease <- "FXS"
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
                      m.gen$TE.random, 
                      m.gen$statistic.random,
                      m.gen$lower.random, 
                      m.gen$upper.random))
  }
}
colnames(output) <- c("Gene", "Pvalue", "logFC", "Statistic", "Lower", "Upper")
outputDF <- as.data.frame(output)
save(outputDF, file = "6. DiseaseMeta/Gene/outputDF_FXS.RData")



################################################################################

# Make gene logFC plot: repeat for DS, DMD, RTT, and FX

################################################################################

disease <- "DMD"

# Load data
load(paste0("6. DiseaseMeta/Gene/outputDF_",disease,".RData"))

# Prepare data
outputDF <- outputDF[!is.na(outputDF$Pvalue),]
outputDF$Gene <- as.character(outputDF$Gene)
geneInfo$GeneID <- as.character(geneInfo$GeneID)
outputDF <- left_join(outputDF, geneInfo, by = c("Gene" = "GeneID"))

# Export statistics in .csv format
exportCSV <- unique(outputDF[,c("Symbol", "Gene","EnsemblGeneID", "Statistic", 
                                "Lower", "Upper","Pvalue")])
colnames(exportCSV) <- c("HGNC", "ENTREZ", "ENSEMBL", "Statistic", "Lower",
                         "Upper", "Pvalue")
write.csv(exportCSV,file = paste0("6. DiseaseMeta/Gene/MetaAnalysis_", disease, ".csv"),
          row.names = FALSE)

# For DS: exclude chromosome 21 genes
#outputDF <- outputDF[outputDF$ChrName != "21",]

# Sort based on P values
outputDF <- arrange(outputDF, by = Pvalue)

# Select gene
selGene <- "5217" #PFN2: DMD
selGene <- "7570" # ZNF22; DS
selGene <- "55133" # SRBD1: FXS
selGene <- "3691" # ITGB4: RTT

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

# Make plot (Vertical)
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


# Make plot (Horizontal)
p <- ggplot(plotDF[plotDF$Disease != "Other",]) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = -BG, ymax = BG), 
            fill = colors[1], alpha = 0.005) +
  geom_hline(yintercept = 0, size = 1, color = "#D9D9D9") +
  geom_segment(aes(x = StudyID, xend = StudyID, y = Lower, yend = Upper, 
                   color = Disease1)) +
  geom_point(aes(x = StudyID, y = logFC, color = Disease1)) +
  facet_grid(cols = vars(`Disease abbreviation`), space = "free", scale = "free") +
  scale_color_manual(values = colors) +
  ylab(expression(log[2]~"FC")) +
  xlab("Datasets") +
  ggtitle(NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

# Save plot
ggsave(p, file = paste0("6. DiseaseMeta/Gene/", disease, "horizontal.png"), 
       height = 2, width = 6.5)


# Export figure source data
sourceData <- plotDF[plotDF$Disease != "Other", c("StudyID", "Disease abbreviation","logFC", "Upper", "Lower")]
colnames(sourceData) <- c("Dataset", "Disease", "logFC", "Upper", "Lower")
write.csv(sourceData, file = "6. DiseaseMeta/SourceData_Figure3B_DMD.csv",
          row.names = FALSE, quote = FALSE)


################################################################################

# Make Manhattan plot: repeat for DS, DMD, RTT, and FX

################################################################################

disease <- "DMD"

# Load data
load(paste0("6. DiseaseMeta/Gene/outputDF_",disease,".RData"))

# Prepare data
outputDF <- outputDF[!is.na(outputDF$Pvalue),]
outputDF$Gene <- as.character(outputDF$Gene)
outputDF$Pvalue <- as.numeric(outputDF$Pvalue)
outputDF$Statistic <- as.numeric(outputDF$Statistic)
geneInfo$GeneID <- as.character(geneInfo$GeneID)
outputDF <- left_join(outputDF, geneInfo, by = c("Gene" = "GeneID"))

outputDF$ChrStart <- as.numeric(outputDF$ChrStart)
outputDF$ChrName[outputDF$ChrAcc == "NC_012920.1"] <- "MT"
outputDF <- outputDF[!is.na(outputDF$ChrName),]
outputDF <- outputDF[!(outputDF$ChrName) %in% c("Y", "MT"),]

outputDF$ChrName <- factor(outputDF$ChrName,
                           levels = c("1", "2", "3", "4", "5","6", "7", "8",
                                      "9", "10", "11", "12", "13", "14", "15",
                                      "16", "17", "18", "19", "20", "21", "22",
                                      "X", "Y", "MT"))
outputDF$Sig <- ifelse(outputDF$Pvalue < 0.05, "Yes", "No")
outputDF$Color <- paste0(outputDF$ChrName, "_", outputDF$Sig)

if (disease == "DMD"){
  mainColor <- "#737373"
  mainColor2 <- "#7570B3"
  sideColor <- "#D9D9D9"
  sideColor2 <- "#DADAEB"
}
if (disease == "DS"){
  mainColor <- "#E6AB02"
  mainColor2 <- "#737373"
  sideColor <- "#feedbd"
  sideColor2 <- "#D9D9D9"
}
if (disease == "FXS"){
  mainColor <- "#737373"
  mainColor2 <- "#E7298A"
  sideColor <- "#D9D9D9"
  sideColor2 <- "#f9cce3"
}
if (disease == "RTT"){
  mainColor <- "#D95F02"
  mainColor2 <- "#737373"
  sideColor <- "#FEE6CE"
  sideColor2 <- "#D9D9D9"
}


colors <- setNames(c(rep(c(mainColor,mainColor2),12),mainColor,
                   rep(c(sideColor,sideColor2),12),sideColor),
                   c("1_Yes", "2_Yes", "3_Yes", "4_Yes", "5_Yes",
                     "6_Yes", "7_Yes", "8_Yes",
                     "9_Yes", "10_Yes", "11_Yes", "12_Yes", 
                     "13_Yes", "14_Yes", "15_Yes",
                     "16_Yes", "17_Yes", "18_Yes", "19_Yes", 
                     "20_Yes", "21_Yes", "22_Yes",
                     "X_Yes", "Y_Yes", "MT_Yes",
                     "1_No", "2_No", "3_No", "4_No", "5_No",
                     "6_No", "7_No", "8_No",
                     "9_No", "10_No", "11_No", "12_No", 
                     "13_No", "14_No", "15_No",
                     "16_No", "17_No", "18_No", "19_No", 
                     "20_No", "21_No", "22_No",
                     "X_No", "Y_No", "MT_No"))

p <- ggplot(outputDF) +
  geom_point(aes(x = as.numeric(ChrStart), y = Statistic, color = Color), 
             size = 1) +
  facet_grid(cols = vars(ChrName), scale = "free") +
  geom_hline(yintercept = 0) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = c(-8,-6,-4,-2,0,2,4,6,8)) +
  xlab(NULL) +
  ylab("t-statistic") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


ggsave(p, file = paste0("6. DiseaseMeta/Gene/", disease, "Manhattan.png"), 
       width = 7, height = 2)



# Export figure source data
sourceData <- outputDF[,c("Gene", "Symbol", "Statistic", "ChrName", "ChrStart")]
colnames(sourceData) <- c("EntrezID", "GeneSymbol", "tstatistic", "Chromosome", "Position")
write.csv(sourceData, file = "6. DiseaseMeta/SourceData_Figure3A_DMD.csv",
          row.names = FALSE, quote = FALSE)





  