# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(clusterProfiler)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("Data/CleanData/statistics_matrix.RData")
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("Data/CleanData/geneInfo.RData")
geneInfo$GeneID <- as.character(geneInfo$GeneID)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Select GO term of interest
selTerm <- "GO:0021694" # Cerebellar Purkinje cell layer formation
selName <- firstup(NES$Description[NES$ID == selTerm])
selGenes <- GOgenes[[selTerm]]

# Select phenotype of interest
selPheno <- "Seizure"
selSamples <- metaData_all$ID[metaData_all[,selPheno] == 1]

# Update matrix
logFC_matrix1 <- logFC_matrix
logFC_matrix1[is.na(logFC_matrix1)] <- 0
SE_matrix1 <- SE_matrix
SE_matrix1[is.na(SE_matrix1)] <- 0

# Calculate OR and pvalues
OR <- rep(NA, length(selGenes))
lw <- rep(NA, length(selGenes))
up <- rep(NA, length(selGenes))
pvalue <- rep(NA, length(selGenes))
for (i in 1:length(selGenes)){
  sel <- selGenes[i]
  name <- geneInfo$Symbol[geneInfo$GeneID == sel]
  plotDF <- data.frame(logFC = logFC_matrix1[which(rownames(logFC_matrix1) == sel),],
                       Lower = logFC_matrix1[which(rownames(logFC_matrix1) == sel),] - 1.96*SE_matrix1[which(rownames(logFC_matrix1) == sel),],
                       Upper = logFC_matrix1[which(rownames(logFC_matrix1) == sel),] + 1.96*SE_matrix1[which(rownames(logFC_matrix1) == sel),],
                       Pheno = factor(ifelse(colnames(logFC_matrix1) %in% selSamples, selPheno, paste0("No ", selPheno)),
                                      levels = c(selPheno, paste0("No ", selPheno))),
                       Dataset = colnames(logFC_matrix1))
  
  plotDF$Sig <- ifelse((plotDF$Lower > 0) | (plotDF$Upper < 0), "Yes", "No")
  
  t <- table(factor(plotDF$Pheno, levels = rev(levels(plotDF$Pheno))), 
             factor(plotDF$Sig, levels = c("No", "Yes")))
  
  test <- fisher.test(t) 
  OR[i] <- test$estimate
  lw[i] <- test$conf.int[1]
  up[i] <- test$conf.int[2]
  pvalue[i] <- test$p.value
  
}
  

# Create dataset for plotting
plotDF <- data.frame(
  OR=log2(OR),
  lw=lw,
  up=up,
  gene = selGenes
)

# Add gene information
plotDF <- inner_join(plotDF, geneInfo[,c(1,2)],
                     by = c("gene"= "GeneID"))

# Set order of genes
plotDF$Symbol <- factor(plotDF$Symbol,
                      levels = c("ATP7A", "GBA1", "SKOR2", "RORA", "LDB1", "HERC1",
                                 "FAIM2", "AGTPBP1", "TTLL1", "WHRN", "CEND1", "LHX5",
                                 "TTC21B", "SLC25A46", "LHX1"))

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 0
to_add <- data.frame(matrix(NA, empty_bar, ncol(plotDF)))
colnames(to_add) <- colnames(plotDF)
data <- rbind(plotDF, to_add)
data <- data %>% arrange(Symbol)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  summarise(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=abs(OR), fill=OR)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=abs(OR), fill=OR), stat="identity", alpha=0.5) +
  geom_bar(aes(x=as.factor(id), y=abs(OR), fill=OR), stat="identity", alpha=0.5) +
  ylim(-2.5,2.5) +  
  labs(fill = expression(log[2]~"OR")) +
  theme_minimal() +
  theme(
    #legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  scale_fill_gradient2(low = "#2171B5", 
                       mid = "white", 
                       high = "#CB181D", 
                       midpoint = 0,
                       limits = c(-1.5,1.5),
                       oob = scales::squish) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=abs(OR)+0.2, label=Symbol, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE )

p

ggsave(p, file = "7. PhenoMeta/GO_Pheno/CircBar_Purkinje.png", width = 7, height = 5)

