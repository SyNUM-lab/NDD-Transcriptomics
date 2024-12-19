# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- str_replace(str_replace(str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}

# Load data
load("4. GSEA/GSEA_GO/GSEAresults_GO.RData")
load("Data/CleanData/metaData_all.RData")


# Update meta data:

# Add disease info
selDis <- c("ASD", "FXS", "DMD", "DS", "RTT")
metaData_all$FXS <- ifelse(metaData_all$`Disease abbreviation` == "FXS",1,0)
metaData_all$DMD <- ifelse(metaData_all$`Disease abbreviation` == "DMD",1,0)
metaData_all$DS <- ifelse(metaData_all$`Disease abbreviation` == "DS",1,0)
metaData_all$RTT <- ifelse(metaData_all$`Disease abbreviation` == "RTT",1,0)


# Add tissue info
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


# Find associations between GSEA significance and phenotype/disease/tissue
pvalues[is.na(pvalues)] <- 1
NES[is.na(NES)] <- 0
pheno <- as.data.frame(metaData_all[,13:26])
sig_p <- matrix(NA,nrow(pvalues), ncol(pheno))
nes_p <- matrix(NA,nrow(pvalues), ncol(pheno))
for (i in 1:nrow(pvalues)) {
  for (c in 1:ncol(pheno)){
    
    sig_p[i,c] <- wilcox.test(-log10(as.numeric(pvalues[i,3:153]))[pheno[,c] == 1], 
                              -log10(as.numeric(pvalues[i,3:153]))[pheno[,c] == 0])$p.value
    
    nes_p[i,c] <- wilcox.test((sign(as.numeric(NES[i,3:153])) * -log10(as.numeric(pvalues[i,3:153])))[pheno[,c] == 1], 
                              (sign(as.numeric(NES[i,3:153])) * -log10(as.numeric(pvalues[i,3:153])))[pheno[,c] == 0])$p.value
  }
}

rownames(sig_p) <- pvalues$ID
rownames(nes_p) <- NES$ID
colnames(sig_p) <- colnames(pheno)
colnames(nes_p) <- colnames(pheno)
save(sig_p, nes_p, file = "7. PhenoMeta/GO_Pheno/GSEAassociations.RData")

# Load (if not already loaded)
load("7. PhenoMeta/GO_Pheno/GSEAassociations.RData")

# Get adj. P value (unsigned -log10 P value)
sig_p_adj <- apply(sig_p,2,function(x) p.adjust(x, method = "fdr"))

# Only look at GO terms that are significant in five datasets
sigGOs <- pvalues$ID[rowSums(pvalues[,3:153]<0.05) > 5]
sig_p_adj <- sig_p_adj[rownames(sig_p_adj) %in% sigGOs,1:8]

# Get adj. P value (signed -log10 P value)
nes_p_adj <- apply(nes_p,2, function(x) p.adjust(x, method = "fdr"))

# Only look at GO terms that are significant in five datasets
sigGOs <- pvalues$ID[rowSums(pvalues[,3:153]<0.05) > 5]
nes_p_adj <- nes_p_adj[rownames(nes_p_adj) %in% sigGOs,1:8]

# Selected term: "Ceberellar Purkinje cell layer formation"
selTerm <- "GO:0021694"
selName <- firstup(NES$Description[NES$ID == selTerm])
selPheno <- "Seizure"

# Prepare data for plotting
plotDF <- data.frame(NES = as.numeric(NES[NES$ID == selTerm,3:153]),
                     Pvalue = as.numeric(pvalues[pvalues$ID == selTerm,3:153]),
                     Pheno = factor(ifelse(pheno[,selPheno] == 1, selPheno, paste0("No ", selPheno)),
                                    levels = c(selPheno, paste0("No ", selPheno))),
                     Dataset = colnames(NES[,3:153]))

# Make plot
p <- ggplot(plotDF) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_point(aes(x = Dataset, y = -log10(Pvalue)*sign(NES), 
                 color = -log10(Pvalue)*sign(NES))) +
  facet_grid(cols = vars(Pheno), scale = "free", space = "free") +
  scale_color_gradient2(low = "#000072", mid = "#FEE6CE", high = "red", midpoint = 0, 
                        trans = "pseudo_log") +
  ylab(expression("Signed -" ~ log[10]~"P value")) +
  xlab("Datasets") +
  ggtitle(selName) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        strip.text = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 12))

# Save plot
ggsave(p, file = paste0("7. PhenoMeta/GO_Pheno/", selName, ".png"), width = 6, height = 4)


################################################################################

# Get genes associated with Cerebellar Purkinje cell layer formation

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/Analysis")

# Load data
load("4. GSEA/GOgenes_BP_ENTREZID_Hs.RData")

# Cerebellar Purkinje cell layer formation
selTerm <- "GO:0021694"
selGenes <- GOgenes[[selTerm]]
write.table(geneInfo$Symbol[geneInfo$GeneID %in% selGenes],
            file = "7. PhenoMeta/GO_Pheno/PurkinjeGenes.txt", quote = FALSE,
            sep = "\t", col.names = FALSE, row.names = FALSE)

# This list of genes can be used as input for GeneMANIA (Cytoscape)