
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(GEOquery)
library(tidyverse)

# Set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")


# https://wires.onlinelibrary.wiley.com/doi/full/10.1002/wcs.1398
NDD_names <- c("Lesch-Nyhan syndrome",
               "Lowe syndrome",
               "Rubinstein-Taybi syndrome",
               "Cornelia de Lange syndrome",
               "Cri du chat syndrome",
               "Galactosaemia",
               "Angelman syndrome",
               "Williams syndrome",
               "Marfan syndrome",
               "Prader-Willi syndrome",
               "Rett syndrome",
               "Phenylketonuria",
               "Duchenne muscular dystrophy",
               "Tuberous sclerosis",
               "Trisomy 18",
               "Velocardiofacial syndrome",
               "Neurofibromatosis type 1",
               "Turner syndrome",
               "XYY",
               "XXX",
               "Noonan syndrome",
               "Fragile X syndrome",
              "Klinefelter syndrome",
               "Fetal alcohol syndrome",
               "Cerebral palsy",
               "Down syndrome",
               "Tourette syndrome",
               "Autism spectrum disorder",
               "Developmental dyscalculia",
               "Attention deficit hyperactivity disorder",
               "Intellectual disability",
               "Developmental dyslexia",
               "Developmental coordination disorder",
               "Specific language impairment",
               "Speech sound disorder")

NDD_names <- tolower(NDD_names)
search <- paste(NDD_names, collapse = '"[All Fields] OR "')
search_final <- paste0('("',search,'") AND "GEO2R"[All Fields] AND "Homo sapiens"[porgn] AND "Expression profiling by high throughput sequencing"[Filter]')

# Export GEO search terms
write.table(search_final,file = "1. PrepareData/SearchQuery.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# [Use search terms in GEO website...]

# Get all study IDs
all_IDs <- read.delim("gds_result.txt", header = FALSE)
all_IDs <- str_replace(all_IDs$V1, "200", "GSE")
all_IDs <- str_replace(all_IDs, "GSE0", "GSE")

# Set working directory to datasets
setwd("D:/RTTproject/GEOData/Data/Datasets")

# Retrieve count and meta data
for (i in 1:length(all_IDs)){
  
  # Get GEO ID
  ID <- all_IDs[i]
  
  # load counts table from GEO
  urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
  path <- paste(urld, paste0("acc=", ID), paste0("file=", ID, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"), sep="&");
  tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
  
  # Only if there are more samples, we want to save the data
  if (ncol(tbl) >= 6){
    
    # Save expression table
    dir.create(ID)
    write.table(tbl, file = paste0(ID,"/",ID, "_expr.tsv"), sep = "\t", quote = FALSE)
    
    # Save meta data
    geoObjectject <- getGEO(ID)
    metaData <- geoObjectject[[1]]@phenoData@data
    metaData$SampleID <- rownames(metaData)
    write.table(metaData, file = paste0(ID,"/",ID, "_meta.tsv"), sep = "\t", quote = FALSE)
  }
  
}

