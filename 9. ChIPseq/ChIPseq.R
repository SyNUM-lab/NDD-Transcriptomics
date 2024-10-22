
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(rtracklayer)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm9)
library(tidyverse)

# Load biomaRt dataset
mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl',
               host="may2009.archive.ensembl.org")

# Load genome annotations for each chromosome
egs = getBM(attributes = c('ensembl_gene_id','external_gene_id',
                           'chromosome_name','start_position',
                           'end_position','strand'), 
            filters='chromosome_name',
            values=1:19,
            mart=mart)

# Get transcription start site
egs$TSS = ifelse(egs$strand == "1", egs$start_position, egs$end_position)

# Save genome annotations
save(egs, file = "9. ChIPseq/egs.RData")

# Read raw data (GSE70957)
D_treat <- rtracklayer::import.wig("Data/ChIPseq/GSE70957_D1_D2merged_for_Wig_treat_afterfiting_all.wig.gz")
D_treat_fil <- D_treat[!(seqnames(D_treat) %in% c("chrM", "chrY", "chrX"))]

# Save data as R object
save(D_treat_fil, file = "9. ChIPseq/D_treat_fil.RData")

#==============================================================================#
# 1) Permutation
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("D:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(rtracklayer)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm9)
library(tidyverse)

# Get sequence info
genome = BSgenome.Mmusculus.UCSC.mm9
si = seqinfo(genome)
si = si[ paste0('chr', c(1:19))]

# Load data
load("9. ChIPseq/D_treat_fil.RData")
load("9. ChIPseq/egs.RData")
load("9. ChIPseq/GOgenes_BP_ENSEMBL_Mm.RData")
load("9. ChIPseq/GOannotation.RData")

# Filter for common genes
for (i in 1:length(GOgenes)){
  GOgenes[[i]] <- intersect(GOgenes[[i]], egs$ensembl_gene_id)
}

# Select GO terms of interest
selTerms <- c("GO:0051402", "GO:0022008")
selNames <- c("Neuron apoptotic process",
              "Neurogenesis")


#==============================================================================#
# Permutation analysis
#==============================================================================#

# Settings:
range <- 3000 # Range around TSS
lo <- 60      # Length out (bin size)
nGenes <- min(unlist(lapply(GOgenes[selTerms], length))) # Number of selected genes
nPerm <- 1000  # Number of permutations
all_genes <- intersect(unique(unlist(GOgenes)), egs$ensembl_gene_id)

permResults <- matrix(NA, nrow = nPerm, ncol = lo)
set.seed(123)
for (p in 858:nPerm){
  
  selGenes <- all_genes[sample(1:length(all_genes),nGenes)]
  egs_sel <- egs[egs$ensembl_gene_id %in% selGenes,]
  
  # Set tiles around TSS of selected genes
  tiles = sapply(1:nrow(egs_sel), function(i)
    if(egs_sel$strand[i] == "1" )
      egs_sel$TSS[i] + seq(-1*range, range-100, length.out=lo)
    else
      egs_sel$TSS[i] + seq(range-100, -1*range, length.out=lo))
  
  tiles = GRanges(tilename = paste(rep(egs_sel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                  seqnames = Rle(rep(paste0('chr', egs_sel$chromosome_name), each=lo)), 
                  ranges = IRanges(start = as.vector(tiles),
                                   width = 100),
                  strand = Rle(rep("*", length(as.vector(tiles)))),
                  seqinfo=si)
  
  # Get bin counts (treat)
  Loci = countOverlaps(tiles, D_treat_fil)
  
  Loci_matrix_sel = matrix(Loci, nrow=nrow(egs_sel), 
                           ncol=lo, byrow=TRUE)
  
  permResults[p,] <- colMeans(Loci_matrix_sel)
  
}
save(permResults, file = "9. ChIPseq/permResults.RData")


#==============================================================================#

#==============================================================================#


# Get mean Lhx1 binding
plotDF_sel <- NULL
for (t in 1:length(selTerms)){
  selGenes <- GOgenes[[selTerms[t]]]
  
  egs_sel <- egs[egs$ensembl_gene_id %in% selGenes,]
  
  # Settings:
  range <- 3000 # Range around TSS
  lo <- 60      # length out (bin size)
  
  
  # Set tiles around TSS of selected genes
  tiles = sapply(1:nrow(egs_sel), function(i)
    if(egs_sel$strand[i] == "1" )
      egs_sel$TSS[i] + seq(-1*range, range-100, length.out=lo)
    else
      egs_sel$TSS[i] + seq(range-100, -1*range, length.out=lo))
  
  tiles = GRanges(tilename = paste(rep(egs_sel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                  seqnames = Rle(rep(paste0('chr', egs_sel$chromosome_name), each=lo)), 
                  ranges = IRanges(start = as.vector(tiles),
                                   width = 100),
                  strand = Rle(rep("*", length(as.vector(tiles)))),
                  seqinfo=si)
  
  # Get bin counts
  Loci = countOverlaps(tiles, D_treat_fil)
  
  # Calculate the mean count
  Loci_matrix_sel = matrix(Loci, nrow=nrow(egs_sel), 
                           ncol=lo, byrow=TRUE)
  selResults_treat <- colMeans(Loci_matrix_sel)
  
  # Combine information into data frame
  temp <- data.frame(Location = seq(-1*range, range, length.out=lo),
                     Treat = selResults_treat,
                     Term = selTerms[t],
                     Name = selNames[t])
  
  plotDF_sel <- rbind.data.frame(plotDF_sel, temp)
}




# Set names of GO terms
plotDF_sel$Name <- paste0(plotDF_sel$Name,"\n(", plotDF_sel$Term, ")")

# Prepare permutation results for plotting
load("9. ChIPseq/permResults.RData")
plotDF <- gather(as.data.frame(permResults))
plotDF$Perm <- rep(1:nPerm,ncol(permResults))
plotDF$key <- rep(seq(-1*range, range, length.out=lo),each = nrow(permResults))
colnames(plotDF) <- c("Location", "Treat", "Perm")

# Set colors
colors <- setNames(c("#A50F15", "#08519C"),
                   c("Neuron apoptotic process\n(GO:0051402)", 
                     "Neurogenesis\n(GO:0022008)"))

# Make plot
p <- ggplot() +
   geom_line(data = plotDF, aes(x = Location, y = Treat, group = Perm), 
             color = "grey", alpha = 0.1) +
  geom_line(data = plotDF_sel, aes(x = Location, y = Treat, color = Name),
            alpha = 1, linewidth = 1.5) +
  xlab("Distance from TSS") +
  ylab("Mean tag count") +
  theme_bw() +
  scale_color_manual(values = colors) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(p, file = "9. ChIPseq/ChIPseqPlot.png", width = 7, height = 5)

# Estimate P value
sum(rowSums(permResults) > sum(plotDF_sel$Treat[plotDF_sel$Term == "GO:0051402"]))
sum(rowSums(permResults) > sum(plotDF_sel$Treat[plotDF_sel$Term == "GO:0022008"]))

