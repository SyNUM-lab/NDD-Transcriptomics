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
library(rrvgo)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Prepare meta data:

# 1) Add disease info
selDis <- c("ASD", "FXS", "DMD", "DS", "RTT")
metaData_all$ASD <- ifelse(metaData_all$`Disease abbreviation` == "ASD",1,0)
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
metaData_all <- arrange(metaData_all, by = ID)
all(colnames(pvalues[,3:153]) == metaData_all$ID)


# Format p values (GSEA)
p_matrix <- as.matrix(pvalues[,3:153])
p_matrix[is.na(p_matrix)] <- 1
rownames(p_matrix) <- pvalues$ID

# calculate adjusted p-values (GSEA)
adjp_matrix <- matrix(p.adjust(p_matrix, method = "fdr"),
                      nrow = nrow(p_matrix), ncol = ncol(p_matrix))
rownames(adjp_matrix) <- rownames(p_matrix)
colnames(adjp_matrix) <- colnames(p_matrix)

scores <- sort(setNames(rowSums(adjp_matrix < 0.05), rownames(adjp_matrix)))
rownames(pvalues) <- pvalues$ID
rownames(NES) <- NES$ID

cbind(pvalues[rev(names(tail(scores,30))), 1:2],rev(tail(scores,30)))

#==============================================================================#
# Make similarity matrix of GO terms
simMatrix <- calculateSimMatrix(
  rownames(adjp_matrix),
  orgdb="org.Hs.eg.db",
  ont = "BP",
  method = "Resnik"
)

# Select GO term with highest number of significant values
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")

save(reducedTerms, file = "5. overallMeta/GSEA/reducedTerms.RData")
#==============================================================================#

load("5. OverallMeta/GSEA/reducedTerms.RData")
selTerms <- reducedTerms[,c("parent", "parentTerm", "score")]
selTerms <- selTerms %>%
  group_by(parent) %>%
  mutate(score = max(score))
rownames(selTerms) <- NULL
selTerms <- unique(selTerms)

# Plot top 30 terms (with most # of significance)
plotTerms <- head(selTerms$parent,30)
pvalues_plot <- pvalues[plotTerms,]
NES_plot <- NES[plotTerms,]

# Prepare data for plotting
plotDF <- gather(pvalues_plot[,3:153])
plotDF$ID <- rep(pvalues_plot$ID,ncol(pvalues_plot)-2)
plotDF$Name <- rep(pvalues_plot$Description,ncol(pvalues_plot)-2)
plotDF$NES <- gather(NES_plot[,3:153])$value
plotDF$value[is.na(plotDF$value)] <- 1
plotDF$NES[is.na(plotDF$NES)] <- 0

# Cluster the GO terms
values_p <- -log10(pvalues_plot[,3:153]) * sign(NES_plot[,3:153])
values_p[is.na(values_p)] <- 0
values_p[values_p > 5] <- 5
values_p[values_p < -5] <- -5
model <- hclust(dist(values_p), "ward.D2")
terms_ordered <- model$labels[model$order]
clusterDF <- data.frame(ID = names(cutree(model, k = 5)),
                        Cluster = cutree(model, k = 5))
save(terms_ordered, clusterDF, file = "5. overallMeta/GSEA/terms_ordered.RData")

# Cluster the datasets
model <- hclust(dist(t(values_p)), "ward.D2")
sets_ordered <- model$labels[model$order]
save(sets_ordered, file = "5. overallMeta/GSEA/sets_ordered.RData")

# Put the terms and datasets in correct order
plotDF <- inner_join(plotDF, clusterDF, by = c("ID" = "ID"))
plotDF$key <- factor(plotDF$key, levels = sets_ordered)
plotDF$ID <- factor(plotDF$ID, levels = terms_ordered)

plotDF$Name <- firstup(plotDF$Name)
plotDF$Name[nchar(plotDF$Name)>50] <- paste0(substring(plotDF$Name[nchar(plotDF$Name)>50],1,47),"...")

plotDF <- arrange(plotDF, by = ID)
plotDF$Name <- factor(plotDF$Name, levels = unique(plotDF$Name))
plotDF$Cluster[plotDF$Cluster == 1] <- "A"
plotDF$Cluster[plotDF$Cluster == 2] <- "B"
plotDF$Cluster[plotDF$Cluster == 3] <- "A"
plotDF$Cluster[plotDF$Cluster == 4] <- "D"
plotDF$Cluster[plotDF$Cluster == 5] <- "C"

# NOTE: cluster 1 and 3 were combined given their similariy in logFCs and 
# included genes

# Make plot
p <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log10(value)*sign(NES)), height = 0.9) +
  facet_grid(rows = vars(Cluster), scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                        trans = "pseudo_log") +
  ylab(NULL) +
  xlab("Datasets") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

panel_colors <- c("#1B9E77","#D95F02", "#E7298A","#E6AB02")

# convert to grob
gp <- ggplotGrob(p)
for(i in 1:4){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}

# Save plot
ggsave(gp, file = "5. overallMeta/GSEA/GSEAplot_all.png", width = 10, height = 6)


# Get legends:
legendPlot <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log10(value)*sign(NES)), height = 0.9) +
  facet_grid(rows = vars(Cluster), scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                       trans = "pseudo_log") +
  ylab(NULL) +
  xlab("Datasets") +
  theme_bw() +
  theme(legend.position = "right",
        strip.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "5. overallMeta/GSEA/GSEAplot_all_legend.png", width = 4, height = 4)


# Export figure source data
sourceData <- plotDF[,c("key", "ID", "Name", "Cluster", "value", "NES")]
colnames(sourceData) <- c("Dataset", "GOID", "GOname", "Cluster", "pvalue", "NES")
write.csv(sourceData, file = "5. overallMeta/SourceData_Figure2A.csv",
          row.names = FALSE, quote = FALSE)

################################################################################

# Phenotype associations

################################################################################

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
library(rrvgo)

# Load data
load("Data/CleanData/topList.RData")
load("Data/CleanData/metaData_all.RData")
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("Data/CleanData/metaData_all.RData")


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
all(colnames(pvalues[,3:153]) == metaData_all$ID)


#==============================================================================#
# Calculate correlations with phenotype
pvalues[is.na(pvalues)] <- 1
NES[is.na(NES)] <- 0
pheno <- as.data.frame(metaData_all[,15:31])
sig_p <- matrix(NA,nrow(pvalues), ncol(pheno))
nes_p <- matrix(NA,nrow(pvalues), ncol(pheno))
sig_sign <- matrix(NA,nrow(pvalues), ncol(pheno))
nes_sign <- matrix(NA,nrow(pvalues), ncol(pheno))
for (i in 1:nrow(pvalues)) {
  for (c in 1:ncol(pheno)){
    sig_p[i,c] <- wilcox.test(-log10(as.numeric(pvalues[i,3:153]))[pheno[,c] == 1], 
                              -log10(as.numeric(pvalues[i,3:153]))[pheno[,c] == 0])$p.value
    sig_sign[i,c] <- median(-log10(as.numeric(pvalues[i,3:153]))[pheno[,c] == 1]) - median(-log10(as.numeric(pvalues[i,3:153]))[pheno[,c] == 0])
    
    nes_p[i,c] <- wilcox.test((sign(as.numeric(NES[i,3:153])) * -log10(as.numeric(pvalues[i,3:153])))[pheno[,c] == 1], 
                              (sign(as.numeric(NES[i,3:153])) * -log10(as.numeric(pvalues[i,3:153])))[pheno[,c] == 0])$p.value
    nes_sign[i,c] <- median((sign(as.numeric(NES[i,3:153])) * -log10(as.numeric(pvalues[i,3:153])))[pheno[,c] == 1]) - median((sign(as.numeric(NES[i,3:153])) * -log10(as.numeric(pvalues[i,3:153])))[pheno[,c] == 0])
    
  }
}

rownames(sig_p) <- pvalues$ID
rownames(nes_p) <- NES$ID
rownames(sig_sign) <- pvalues$ID
rownames(nes_sign) <- NES$ID
colnames(sig_p) <- colnames(pheno)
colnames(nes_p) <- colnames(pheno)
colnames(sig_sign) <- colnames(pheno)
colnames(nes_sign) <- colnames(pheno)
save(nes_p, sig_p, nes_sign, sig_sign, file = "5. overallMeta/GSEA/GO_pheno_correlations1.RData")
#==============================================================================#

# Load data
load("5. overallMeta/GSEA/GO_pheno_correlations1.RData")
load("5. overallMeta/GSEA/terms_ordered.RData")

# Plot correlations with signed GSEA P value:

# Prepare data for plotting
plotDF <- gather(as.data.frame(-log10(nes_p[terms_ordered,1:17])))
plotDF$sign <- sign(gather(as.data.frame(sign(nes_sign[terms_ordered,1:17])))$value)
plotDF$ID <- rep(terms_ordered,17)
plotDF <- inner_join(plotDF, clusterDF, by = c("ID" = "ID"))
plotDF <- inner_join(plotDF,pvalues[,1:2], by = c("ID" = "ID"))
plotDF$Description <- firstup(plotDF$Description)
plotDF$Description[nchar(plotDF$Description)>50] <- paste0(substring(plotDF$Description[nchar(plotDF$Description)>50],1,47),"...")

# Set order and names
plotDF$ID <- factor(plotDF$ID,levels = terms_ordered)
plotDF <- arrange(plotDF, by = ID)
plotDF$Description <- factor(plotDF$Description, levels = unique(plotDF$Description))
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
plotDF$Cluster[plotDF$Cluster == 1] <- "A"
plotDF$Cluster[plotDF$Cluster == 2] <- "B"
plotDF$Cluster[plotDF$Cluster == 3] <- "A"
plotDF$Cluster[plotDF$Cluster == 4] <- "D"
plotDF$Cluster[plotDF$Cluster == 5] <- "C"
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
  geom_tile(aes(x = key, y = Description, fill = value), height = 0.9) +
  geom_text(data = plotDF[plotDF$value > -log10(0.05),], 
            aes(x = key, y = Description, label = sign), color = "white") +
  scale_fill_gradient(low = "#FCFBFD", high = "#4A1486", limits = c(0,3),
                      trans = "pseudo_log", oob = scales::squish) +
  facet_grid(rows = vars(Cluster), 
             cols = vars(varGroup),
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

panel_colors <- c("#1B9E77","#D95F02", "#E7298A","#E6AB02")

# convert to grob
gp <- ggplotGrob(p)
for(i in 1:4){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}

# Save plot
ggsave(gp, file = "5. overallMeta/GSEA/GSEAplot_all_cor_signed.png", width = 6, height = 8)




# Plot correlations with unsigned GSEA P value:

# Prepare data for plotting
plotDF <- gather(as.data.frame(-log10(sig_p[terms_ordered,1:17])))
plotDF$sign <- sign(gather(as.data.frame(sign(nes_sign[terms_ordered,1:17])))$value)
plotDF$ID <- rep(terms_ordered,17)
plotDF <- inner_join(plotDF, clusterDF, by = c("ID" = "ID"))
plotDF <- inner_join(plotDF,pvalues[,1:2], by = c("ID" = "ID"))
plotDF$Description <- firstup(plotDF$Description)
plotDF$Description[nchar(plotDF$Description)>50] <- paste0(substring(plotDF$Description[nchar(plotDF$Description)>50],1,47),"...")

# Set order and names
plotDF$ID <- factor(plotDF$ID,levels = terms_ordered)
plotDF <- arrange(plotDF, by = ID)
plotDF$Description <- factor(plotDF$Description, levels = unique(plotDF$Description))
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
plotDF$Cluster[plotDF$Cluster == 1] <- "A"
plotDF$Cluster[plotDF$Cluster == 2] <- "B"
plotDF$Cluster[plotDF$Cluster == 3] <- "A"
plotDF$Cluster[plotDF$Cluster == 4] <- "D"
plotDF$Cluster[plotDF$Cluster == 5] <- "C"
plotDF$sign <- ifelse(plotDF$sign == 1, "#", "0")

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
  geom_tile(aes(x = key, y = Description, fill = value), height = 0.9) +
  geom_text(data = plotDF[(plotDF$value > -log10(0.05)) & (plotDF$sign == "#"),], 
            aes(x = key, y = Description, label = sign), color = "white") +
  scale_fill_gradient(low = "#F7FCF5", high = "#238B45", limits = c(0,3),
                      trans = "pseudo_log", oob = scales::squish) +
  facet_grid(rows = vars(Cluster), 
             cols = vars(varGroup),
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

panel_colors <- c("#1B9E77","#D95F02", "#E7298A","#E6AB02")

# convert to grob
gp <- ggplotGrob(p)
for(i in 1:4){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}

# Save plot
ggsave(gp, file = "5. overallMeta/GSEA/GSEAplot_all_cor_unsigned.png", width = 6, height = 8)


