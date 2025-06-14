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
library(meta)
library(rrvgo)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), 
                   "RRNA", "rRNA")
  return(x)
}

################################################################################

# Find relevant GO terms: repeat for DS, DMD, RTT, and FX

################################################################################
simMatrix <- calculateSimMatrix(
  pvalues$ID,
  orgdb="org.Hs.eg.db",
  ont = "BP",
  method = "Resnik"
)
save(simMatrix, file = "6. DiseaseMeta/GO/simMatrix.RData")

#==============================================================================#
# DMD
#==============================================================================#

disease <- "DMD"

# Select terms without DMD
terms <- names(GOgenes)[unlist(lapply(GOgenes, function(x) sum(str_detect(x, "1756"))==0))]

# Select all DMD datasets
samples <- metaData_all$ID[metaData_all$`Disease abbreviation` == disease]
NES_fil <- NES[NES$ID %in% terms,c("ID", "Description", samples)]
NES_fil[is.na(NES_fil)] <- 1
pvalues_fil <- pvalues[pvalues$ID %in% terms,c("ID", "Description", samples)]
pvalues_fil[is.na(pvalues_fil)] <- 1

# number of significant datasets
test_p <- pvalues_fil[,3:ncol(pvalues_fil)] < 0.05

# Number of downregulated datasets
test_neg <- NES_fil[,3:ncol(NES_fil)] < 0

# Number for upregulated datasets
test_pos <- NES_fil[,3:ncol(NES_fil)] > 0

# Combine into dataframe
resultDF <- data.frame(ID = pvalues_fil$ID,
                       Description = pvalues_fil$Description,
                       total = rowSums(test_p),
                       neg = rowSums(test_p & test_neg),
                       pos = rowSums(test_p & test_pos)
)

# Get number of significant datasets with consistent direction of effect
resultDF$value <- apply(resultDF[,4:5],1,max)


# Select GO term with highest number of significant values
reducedTerms <- reduceSimMatrix(simMatrix[rownames(simMatrix) %in% resultDF$ID,
                                          colnames(simMatrix) %in% resultDF$ID],
                                scores = sort(setNames(resultDF$value[resultDF$ID %in% rownames(simMatrix)],
                                                  resultDF$ID[resultDF$ID %in% rownames(simMatrix)])),
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")




#==============================================================================#
# DS
#==============================================================================#

disease <- "DS"

# Select all DMD datasets
samples <- metaData_all$ID[metaData_all$`Disease abbreviation` == disease]
NES_fil <- NES[NES$ID %in% terms,c("ID", "Description", samples)]
NES_fil[is.na(NES_fil)] <- 1
pvalues_fil <- pvalues[pvalues$ID %in% terms,c("ID", "Description", samples)]
pvalues_fil[is.na(pvalues_fil)] <- 1

# number of significant datasets
test_p <- pvalues_fil[,3:ncol(pvalues_fil)] < 0.05

# Number of downregulated datasets
test_neg <- NES_fil[,3:ncol(NES_fil)] < 0

# Number for upregulated datasets
test_pos <- NES_fil[,3:ncol(NES_fil)] > 0

# Combine into dataframe
resultDF <- data.frame(ID = pvalues_fil$ID,
                       Description = pvalues_fil$Description,
                       total = rowSums(test_p),
                       neg = rowSums(test_p & test_neg),
                       pos = rowSums(test_p & test_pos)
)

# Get number of significant datasets with consistent direction of effect
resultDF$value <- apply(resultDF[,4:5],1,max)


# Select GO term with highest number of significant values
reducedTerms <- reduceSimMatrix(simMatrix[rownames(simMatrix) %in% resultDF$ID,
                                          colnames(simMatrix) %in% resultDF$ID],
                                scores = sort(setNames(resultDF$value[resultDF$ID %in% rownames(simMatrix)],
                                                       resultDF$ID[resultDF$ID %in% rownames(simMatrix)])),
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")



#==============================================================================#
# FXS
#==============================================================================#

disease <- "FXS"
terms <- names(GOgenes)[unlist(lapply(GOgenes, function(x) sum(str_detect(x, "2332"))==0))]

# Select all DMD datasets
samples <- metaData_all$ID[metaData_all$`Disease abbreviation` == disease]
NES_fil <- NES[NES$ID %in% terms,c("ID", "Description", samples)]
NES_fil[is.na(NES_fil)] <- 1
pvalues_fil <- pvalues[pvalues$ID %in% terms,c("ID", "Description", samples)]
pvalues_fil[is.na(pvalues_fil)] <- 1

# number of significant datasets
test_p <- pvalues_fil[,3:ncol(pvalues_fil)] < 0.05

# Number of downregulated datasets
test_neg <- NES_fil[,3:ncol(NES_fil)] < 0

# Number for upregulated datasets
test_pos <- NES_fil[,3:ncol(NES_fil)] > 0

# Combine into dataframe
resultDF <- data.frame(ID = pvalues_fil$ID,
                       Description = pvalues_fil$Description,
                       total = rowSums(test_p),
                       neg = rowSums(test_p & test_neg),
                       pos = rowSums(test_p & test_pos)
)

# Get number of significant datasets with consistent direction of effect
resultDF$value <- apply(resultDF[,4:5],1,max)


# Select GO term with highest number of significant values
reducedTerms <- reduceSimMatrix(simMatrix[rownames(simMatrix) %in% resultDF$ID,
                                          colnames(simMatrix) %in% resultDF$ID],
                                scores = sort(setNames(resultDF$value[resultDF$ID %in% rownames(simMatrix)],
                                                       resultDF$ID[resultDF$ID %in% rownames(simMatrix)])),
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")




#==============================================================================#
# RTT
#==============================================================================#

disease <- "RTT"

# Select all DMD datasets
samples <- metaData_all$ID[metaData_all$`Disease abbreviation` == disease]
NES_fil <- NES[NES$ID %in% terms,c("ID", "Description", samples)]
NES_fil[is.na(NES_fil)] <- 1
pvalues_fil <- pvalues[pvalues$ID %in% terms,c("ID", "Description", samples)]
pvalues_fil[is.na(pvalues_fil)] <- 1

# number of significant datasets
test_p <- pvalues_fil[,3:ncol(pvalues_fil)] < 0.05

# Number of downregulated datasets
test_neg <- NES_fil[,3:ncol(NES_fil)] < 0

# Number for upregulated datasets
test_pos <- NES_fil[,3:ncol(NES_fil)] > 0

# Combine into dataframe
resultDF <- data.frame(ID = pvalues_fil$ID,
                       Description = pvalues_fil$Description,
                       total = rowSums(test_p),
                       neg = rowSums(test_p & test_neg),
                       pos = rowSums(test_p & test_pos)
)

# Get number of significant datasets with consistent direction of effect
resultDF$value <- apply(resultDF[,4:5],1,max)


# Select GO term with highest number of significant values
reducedTerms <- reduceSimMatrix(simMatrix[rownames(simMatrix) %in% resultDF$ID,
                                          colnames(simMatrix) %in% resultDF$ID],
                                scores = sort(setNames(resultDF$value[resultDF$ID %in% rownames(simMatrix)],
                                                       resultDF$ID[resultDF$ID %in% rownames(simMatrix)])),
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")



################################################################################

# Make plot: repeat for each GO term

################################################################################

disease <- "DMD"

# Select GO term: 

# DMD
GOterm <- "kidney development"
GOterm <- "cardiac chamber development"
GOterm <- "cell growth"

# DS
GOterm <- "detoxification"
GOterm <- "humoral immune response"

# RTT
GOterm <- "wound healing"
GOterm <- "axonogenesis"

# FXS:
GOterm <- "positive regulation of interleukin-1 beta production"
GOterm <- "cell-cell adhesion via plasma-membrane adhesion molecules"
GOterm <- "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction"

# Make data frame for disease plotting
pvalues[is.na(pvalues)] <- 1
NES[is.na(NES)] <- 0
plotDF <- data.frame(Pvalue = -log10(as.numeric(pvalues[which(pvalues$Description == GOterm)[1], 3:153])),
                     NES = as.numeric(NES[which(NES$Description == GOterm)[1], 3:153]),
                     StudyID = colnames(pvalues)[3:153])
plotDF$signedPvalue <- plotDF$Pvalue*sign(plotDF$NES)

plotDF <- inner_join(plotDF, metaData_all, by = c("StudyID" = "ID"))
plotDF$`Disease abbreviation`[!(plotDF$`Disease abbreviation` %in% c("RTT", "FXS",
                                                                     "DMD", "DS"))] <- "Other"
plotDF$Disease[plotDF$`Disease abbreviation` == "Other"] <- "Other"
plotDF$Disease1 <- ifelse(plotDF$`Disease abbreviation` == disease, "a","b")
plotDF$StudyID <- factor(plotDF$StudyID)

# Set colors
diseaseColors <- setNames(c("#D95F02", "#7570B3", "#E7298A", "#E6AB02"),
                          c("RTT", "DMD", "FXS", "DS"))
colors <- c(diseaseColors[disease], "#BDBDBD")
names(colors) <- c("a", "b")
plotDF$BG <- ifelse(plotDF$`Disease abbreviation`==disease,Inf,0)

# Make plot
p <- ggplot(plotDF[plotDF$Disease != "Other",]) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, 
                xmin = -BG, xmax = BG), 
            fill = colors[1], alpha = 0.005) +
  geom_vline(xintercept = 0, size = 1, color = "#D9D9D9") +
  geom_vline(xintercept = -log10(0.05), size = 0.5, color = "#525252",
             linetype = "dashed") +
  geom_vline(xintercept = log10(0.05), size = 0.5, color = "#525252",
             linetype = "dashed") +
  geom_bar(aes(y = StudyID, x = signedPvalue, fill = Disease1),
           stat = "identity", position = position_dodge()) +
  facet_grid(rows = vars(`Disease abbreviation`), space = "free", scale = "free") +
  scale_fill_manual(values = colors) +
  xlab(expression("Signed -" ~log[10]~"P value")) +
  ylab("Datasets") +
  ggtitle(firstup(GOterm)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5))

# Save plot
ggsave(p, file = paste0("6. DiseaseMeta/GO/", disease, "_GSEA_small3.png"), 
       width = 2, height = 7)


# Export figure source data
sourceData <- plotDF[plotDF$Disease != "Other",c("StudyID", "Disease abbreviation", "Pvalue", "NES")]
colnames(sourceData) <- c("Dataset", "Disease", "Pvalue", "NES")
write.csv(sourceData, file = "6. DiseaseMeta/SourceData_Figure3D_DMD_2.csv",
          row.names = FALSE, quote = FALSE)
