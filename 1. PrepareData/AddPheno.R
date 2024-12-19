
###############################################################################

# Add phenotypes to meta data

##############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)

# Load in phenotype information
phenoFiles <- list.files("Data/Phenotypes")
phenoList <- list()
phenoDF <- NULL

for (p in phenoFiles){
  phenoList[[p]] <- data.table::fread(paste0("Data/Phenotypes/",p))
  
  temp <- data.frame(Phenotype = p,
                        id = phenoList[[p]]$id,
                        name = phenoList[[p]]$name)
  
  phenoDF <- rbind.data.frame(phenoDF, temp)
}

# Select relevant phenotypes
selPheno <- c("seizure", "intellectual disability", "global developmental delay",
  "hypotonia", "scoliosis", "microcephaly", "autis", "gait ataxia")

phenoDF$name <- tolower(phenoDF$name)
phenoMat <- matrix(0,nrow = length(unique(phenoDF$Phenotype)), ncol = length(selPheno))
rownames(phenoMat) <- unique(phenoDF$Phenotype)
colnames(phenoMat) <- selPheno
for (d in 1:nrow(phenoMat)){
  for (p in 1:ncol(phenoMat)){
    if (sum(str_detect(tolower(phenoList[[rownames(phenoMat)[d]]]$name), colnames(phenoMat)[p])) > 0){
      phenoMat[d,p] <- 1
    }
  }
}
colnames(phenoMat) <- c("Seizure", "Intellectual disability", "Global developmental delay",
                        "Hypotonia", "Scoliosis", "Microcephaly", "Autism/Autistic Behavior", "Gait ataxia")
phenoDF_complete <- as.data.frame(phenoMat)
phenoDF_complete$Disease <- rownames(phenoDF_complete)


# Load meta data
metaData <- readxl::read_excel("Data/CleanData/MetaData_clean.xlsx")

# Make sure that each disease has a unique row
diseases <- unique(metaData[,c("Disease","Disease abbreviation","OMIM ID")])
plotDF <- NULL
for (i in 1:nrow(diseases)){
  if (str_detect(diseases$`Disease abbreviation`[i], "; ")){
    temp <- data.frame(
      Disease = str_split(diseases$`Disease`[i], "; ")[[1]],
      `Disease abbreviation` = str_split(diseases$`Disease abbreviation`[i], "; ")[[1]],
      `OMIM ID` = str_split(diseases$`OMIM ID`[i], "; ")[[1]],
      FullDisease = diseases$`Disease`[i]
    )
  } else{
    temp <- data.frame(
      Disease = diseases$`Disease`[i],
      `Disease abbreviation` = diseases$`Disease abbreviation`[i],
      `OMIM ID` = diseases$`OMIM ID`[i],
      FullDisease = diseases$`Disease`[i]
    )
  }
  
  plotDF <- rbind.data.frame(plotDF, temp)
}

# Add phenotype data
plotDF_all <- inner_join(plotDF, phenoDF_complete, by = c("Disease.abbreviation" = "Disease"))
plotDF$`Disease.abbreviation`[which(!(plotDF$`Disease.abbreviation` %in% phenoDF_complete$Disease))]
test <- plotDF_all %>%
  group_by(FullDisease) %>%
  mutate(Seizure = max(Seizure),
         `Intellectual disability` = max(`Intellectual disability`),
         `Global developmental delay` = max(`Global developmental delay`),
         Hypotonia = max(Hypotonia),
         Scoliosis = max(Scoliosis),
         Microcephaly = max(Microcephaly),
         `Autism/Autistic Behavior` = max(`Autism/Autistic Behavior`)
         )
test <- unique(test[,-c(1,2,3)])

# Make final meta data
metaData_all <- inner_join(metaData, test, by = c("Disease" = "FullDisease"))
save(metaData_all, file = "Data/CleanData/metaData_all.RData")



