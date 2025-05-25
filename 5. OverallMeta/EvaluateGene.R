# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

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

pvalue_matrix <- pvalue_matrix[,metaData_all$ID]
logFC_matrix <- logFC_matrix[,metaData_all$ID]
SE_matrix <- SE_matrix[,metaData_all$ID]
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
testDF <- inner_join(testDF, clusterDF, by = c("GO" = "ID"))

# Change cluster names
testDF$Cluster[testDF$Cluster == 1] <- "A"
testDF$Cluster[testDF$Cluster == 2] <- "B"
testDF$Cluster[testDF$Cluster == 3] <- "A"
testDF$Cluster[testDF$Cluster == 4] <- "D"
testDF$Cluster[testDF$Cluster == 5] <- "C"

testDF$gene[testDF$Cluster=="A"][which.max(testDF$value[testDF$Cluster == "A"])]
testDF$gene[testDF$Cluster=="B"][which.max(testDF$value[testDF$Cluster == "B"])]
testDF$gene[testDF$Cluster=="C"][which.max(testDF$value[testDF$Cluster == "C"])]
testDF$gene[testDF$Cluster=="D"][which.max(testDF$value[testDF$Cluster == "D"])]

# Select gene of interest
selGene <- "3106" # cluster A: 45 HLA-B
selGene <- "6451" # cluster B: 36 SH3BGRL
selGene <- "6622" # cluster C 39 SNCA
selGene <- "351" # cluster D: 44 APP

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

# Phenotype associations

################################################################################

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
library(rrvgo)

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

# Prepare meta data:

# 1) Add disease info
selDis <- c("FXS", "DMD", "DS", "RTT")
metaData_all$FXS <- ifelse(metaData_all$`Disease abbreviation` == "FXS",1,0)
metaData_all$DMD <- ifelse(metaData_all$`Disease abbreviation` == "DMD",1,0)
metaData_all$DS <- ifelse(metaData_all$`Disease abbreviation` == "DS",1,0)
metaData_all$RTT <- ifelse(metaData_all$`Disease abbreviation` == "RTT",1,0)

# 2) Add tissue info
neural_tissue <- c("Astrocyte", "BMEC; Astrocyte; Neuron", "Brain organoid",
                   "Cerebral organoid", "Cortical organoid", "Cortical progenitor cell", 
                   "Early neuroectoderm", "iPSC; Neuroepithelium; Neuron", "Late neuroectoderm",
                   "Neuron", "NPC", "Telencephalic organoid", "Brain", "Cerebral granule cell",
                   "Neuroepithelial stem cell"
)
immune_tissue <- c("White blood cell", "PBMC", "Leukocyte", "Lymphocyte", 
                   "Lymphoblastoid cell", "B cell",
                   "Blood", "Microglia cell", "Microglia-like cell")

metaData_all$Neural <- ifelse(metaData_all$Tissue %in% neural_tissue,1,0)
metaData_all$Immune <- ifelse(metaData_all$Tissue %in% immune_tissue,1,0)
metaData_all$Other <- ifelse(metaData_all$Tissue %in% c(neural_tissue, immune_tissue),0,1)

# Add model system info
metaData_all$inVitro<- ifelse(metaData_all$System == "Primary",0,1)
metaData_all$nonInVitro <- ifelse(metaData_all$System == "Primary",1,0)

# sort meta data
metaData_all <- arrange(metaData_all, by = ID)
all(colnames(logFC_matrix) == metaData_all$ID)

#==============================================================================#
# Calculate correlations with phenotype
logFC_matrix_copy <- logFC_matrix
logFC_matrix <- logFC_matrix[rownames(logFC_matrix) %in% c("351", "6622", "6451", "3106"),]
logFC_matrix[is.na(logFC_matrix)] <- 0
pheno <- as.data.frame(metaData_all[,15:31])
sign_p <- matrix(NA,nrow(logFC_matrix), ncol(pheno))
unsign_p <- matrix(NA,nrow(logFC_matrix), ncol(pheno))
sign_dir <- matrix(NA,nrow(logFC_matrix), ncol(pheno))
unsign_dir <- matrix(NA,nrow(logFC_matrix), ncol(pheno))
for (i in 1:nrow(logFC_matrix)) {
  for (c in 1:ncol(pheno)){
    sign_p[i,c] <- wilcox.test(logFC_matrix[i,pheno[,c] == 1], 
                               logFC_matrix[i,pheno[,c] == 0])$p.value
    sign_dir[i,c] <- median(logFC_matrix[i,pheno[,c] == 1] - logFC_matrix[i,pheno[,c] == 0])
    
    unsign_p[i,c] <- wilcox.test(abs(logFC_matrix[i,pheno[,c] == 1]), 
                                 abs(logFC_matrix[i,pheno[,c] == 0]))$p.value
    unsign_dir[i,c] <- median(abs(logFC_matrix[i,pheno[,c] == 1]) - abs(logFC_matrix[i,pheno[,c] == 0]))
    
  }
}

rownames(sign_p) <- rownames(logFC_matrix)
rownames(unsign_p) <- rownames(logFC_matrix)
rownames(sign_dir) <- rownames(logFC_matrix)
rownames(unsign_dir) <- rownames(logFC_matrix)
colnames(sign_p) <- colnames(pheno)
colnames(unsign_p) <- colnames(pheno)
colnames(sign_dir) <- colnames(pheno)
colnames(unsign_dir) <- colnames(pheno)
save(unsign_p, sign_p, unsign_dir, sign_dir, file = "5. overallMeta/GSEA/Gene_pheno_correlations1.RData")
#==============================================================================#

# Load data
load("5. overallMeta/GSEA/Gene_pheno_correlations1.RData")

# Plot correlations with signed GSEA P value:

# Prepare data for plotting
plotDF <- gather(as.data.frame(-log10(sign_p[,1:17])))
plotDF$sign <- sign(gather(as.data.frame(sign(sign_dir[,1:17])))$value)
plotDF$ID <- as.character(rep(rownames(sign_p),17))
plotDF <- inner_join(plotDF, unique(geneInfo[,1:2]), by = c("ID" = "GeneID"))


# Set order and names
plotDF$key <- apply(plotDF, 1, FUN = function(x) {switch(x[1],
                                                         "ASD" = "Autism Spectrum Disorder",
                                                         "DMD" = "Duchenne Muscular Dystrophy",
                                                         "DS" = "Down Syndrome",
                                                         "FXS" = "Fragile X Syndrome",
                                                         "RTT" = "Rett Sydrome",
                                                         "Autism/Autistic Behavior" = "Autistic Behavior",
                                                         "Hypotonia" = "Hypotonia",
                                                         "Gait ataxia" = "Gait Ataxia",
                                                         "Seizure" = "Seizure",
                                                         "Intellectual disability" = "Intellectual Disability",
                                                         "Global developmental delay" = "Global Developmental Delay",
                                                         "Scoliosis" = "Scoliosis",
                                                         "Microcephaly" = "Microcephaly",
                                                         "Neural" = "Neural",
                                                         "Immune" = "Immune",
                                                         "Other" = "Other",
                                                         "inVitro" = "In vitro",
                                                         "nonInVitro" = "Non-in vitro"
)})
plotDF$key <- factor(plotDF$key, 
                     levels = c(             "Autism Spectrum Disorder",
                                             "Duchenne Muscular Dystrophy",
                                             "Down Syndrome",
                                             "Fragile X Syndrome",
                                             "Rett Sydrome",
                                             "Autistic Behavior",
                                             "Hypotonia",
                                             "Gait Ataxia",
                                             "Seizure",
                                             "Intellectual Disability",
                                             "Global Developmental Delay",
                                             "Scoliosis",
                                             "Microcephaly",
                                             "Neural",
                                             "Immune",
                                             "Other",
                                             "In vitro",
                                             "Non-in vitro"))

plotDF$adjp <- p.adjust(10^((plotDF$value)*-1), method = "fdr")
plotDF$sign <- ifelse(plotDF$sign == -1, "-", "+")

plotDF$varGroup <- "Disease"
plotDF$varGroup[plotDF$key %in% c("Autistic Behavior",
                                  "Hypotonia",
                                  "Gait Ataxia",
                                  "Seizure",
                                  "Intellectual Disability",
                                  "Global Developmental Delay",
                                  "Scoliosis",
                                  "Microcephaly")] <- "Phenotype"
plotDF$varGroup[plotDF$key %in% c("Neural",
                                  "Immune",
                                  "Other")] <- "Tissue"
plotDF$varGroup[plotDF$key %in% c("In vitro", "Non-in vitro")] <- "System"

# Make plot
p <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Symbol, fill = value), height = 0.9) +
  geom_text(data = plotDF[plotDF$value > -log10(0.05),], 
            aes(x = key, y = Symbol, label = sign), color = "white") +
  scale_fill_gradient(low = "#FCFBFD", high = "#4A1486", limits = c(0,5),
                      trans = "pseudo_log", oob = scales::squish) +
  facet_grid(cols = vars(varGroup),
             space = "free", scale = "free") +
  ggtitle("Signed") +
  ylab(NULL) +
  xlab(NULL) +
  labs(fill = expression(-log[10] ~ "P value  ")) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_blank(),
        strip.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = "5. overallMeta/GSEA/Gene_cor_signed.png", width = 4, height = 4)




# Load data
load("5. overallMeta/GSEA/Gene_pheno_correlations1.RData")

# Plot correlations with signed GSEA P value:

# Prepare data for plotting
plotDF <- gather(as.data.frame(-log10(unsign_p[,1:17])))
plotDF$sign <- sign(gather(as.data.frame(sign(unsign_dir[,1:17])))$value)
plotDF$ID <- as.character(rep(rownames(unsign_p),17))
plotDF <- inner_join(plotDF, unique(geneInfo[,1:2]), by = c("ID" = "GeneID"))


# Set order and names
plotDF$key <- apply(plotDF, 1, FUN = function(x) {switch(x[1],
                                                         "ASD" = "Autism Spectrum Disorder",
                                                         "DMD" = "Duchenne Muscular Dystrophy",
                                                         "DS" = "Down Syndrome",
                                                         "FXS" = "Fragile X Syndrome",
                                                         "RTT" = "Rett Sydrome",
                                                         "Autism/Autistic Behavior" = "Autistic Behavior",
                                                         "Hypotonia" = "Hypotonia",
                                                         "Gait ataxia" = "Gait Ataxia",
                                                         "Seizure" = "Seizure",
                                                         "Intellectual disability" = "Intellectual Disability",
                                                         "Global developmental delay" = "Global Developmental Delay",
                                                         "Scoliosis" = "Scoliosis",
                                                         "Microcephaly" = "Microcephaly",
                                                         "Neural" = "Neural",
                                                         "Immune" = "Immune",
                                                         "Other" = "Other",
                                                         "inVitro" = "In vitro",
                                                         "nonInVitro" = "Non-in vitro"
)})
plotDF$key <- factor(plotDF$key, 
                     levels = c(             "Autism Spectrum Disorder",
                                             "Duchenne Muscular Dystrophy",
                                             "Down Syndrome",
                                             "Fragile X Syndrome",
                                             "Rett Sydrome",
                                             "Autistic Behavior",
                                             "Hypotonia",
                                             "Gait Ataxia",
                                             "Seizure",
                                             "Intellectual Disability",
                                             "Global Developmental Delay",
                                             "Scoliosis",
                                             "Microcephaly",
                                             "Neural",
                                             "Immune",
                                             "Other",
                                             "In vitro",
                                             "Non-in vitro"))

plotDF$adjp <- p.adjust(10^((plotDF$value)*-1), method = "fdr")
plotDF$sign <- ifelse(plotDF$sign == 1, "#", " ")

plotDF$varGroup <- "Disease"
plotDF$varGroup[plotDF$key %in% c("Autistic Behavior",
                                  "Hypotonia",
                                  "Gait Ataxia",
                                  "Seizure",
                                  "Intellectual Disability",
                                  "Global Developmental Delay",
                                  "Scoliosis",
                                  "Microcephaly")] <- "Phenotype"
plotDF$varGroup[plotDF$key %in% c("Neural",
                                  "Immune",
                                  "Other")] <- "Tissue"
plotDF$varGroup[plotDF$key %in% c("In vitro", "Non-in vitro")] <- "System"

# Make plot
p <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Symbol, fill = value), height = 0.9) +
  geom_text(data = plotDF[plotDF$value > -log10(0.05),], 
            aes(x = key, y = Symbol, label = sign), color = "white") +
  scale_fill_gradient(low = "#F7FCF5", high = "#238B45", limits = c(0,5),
                      trans = "pseudo_log", oob = scales::squish) +
  facet_grid(cols = vars(varGroup),
             space = "free", scale = "free") +
  ggtitle("Unsigned") +
  ylab(NULL) +
  xlab(NULL) +
  labs(fill = expression(-log[10] ~ "P value  ")) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_blank(),
        strip.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = "5. overallMeta/GSEA/Gene_cor_unsigned.png", width = 4, height = 4)




################################################################################

# Network

################################################################################

# Link genes to GO term
networkDF <- NULL
netGenes <- c("351", "6622", "3106", "6451")
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
  title = "GO network2",
  collection = "Genes"
)

