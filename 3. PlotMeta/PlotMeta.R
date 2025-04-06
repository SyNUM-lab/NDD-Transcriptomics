

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# set working directory
setwd("E:/RTTproject/GEOData/NDD-Transcriptomics")

# Load packages
library(tidyverse)
library(ggridges)
library(ggbreak)

# Load data
metaData <- readxl::read_excel("Data/CleanData/MetaData_clean.xlsx")

################################################################################

# barchart of samples sizes

################################################################################

# Unique studies/datasets
nStudies <- unique(str_remove(metaData$ID, "-.*"))

# Prepare data for plotting
plotData <- data.frame(n = c(metaData$`N - disease`, metaData$`N - control`),
                       group = c(rep("Cases (n = 1362)", 151), rep("Controls (n = 1026)", 151)))

# Set colors
colors <- setNames(c("#CB181D", "#2171B5"),
                   c("Cases (n = 1362)", "Controls (n = 1026)"))

colors <- setNames(c("#54278F", "#9E9AC8"),
                   c("Cases (n = 1362)", "Controls (n = 1026)"))

bin_lower <- c(3:9,10,15,20)
bin_upper <- c(3:9,14,19,Inf)
test <- data.frame(table(plotData$n[plotData$group == "Cases (n = 1362)"]))
test$Var1 <- as.numeric(as.character(test$Var1))
binCount_case <- rep(NA, length(bin_lower))
for (i in 1:length(bin_lower)){
  binCount_case[i] <- sum(test$Freq[(test$Var1 >= bin_lower[i]) & (test$Var1 <= bin_upper[i])])
  
}

test <- data.frame(table(plotData$n[plotData$group == "Controls (n = 1026)"]))
test$Var1 <- as.numeric(as.character(test$Var1))
binCount_control <- rep(NA, length(bin_lower))
for (i in 1:length(bin_lower)){
  binCount_control[i] <- sum(test$Freq[(test$Var1 >= bin_lower[i]) & (test$Var1 <= bin_upper[i])])
  
}

plotDF <- data.frame(
  Name = factor(rep(c("3", "4", "5", "6", "7", "8", "9", "10-14", "15-19", "\u2265 20"),2),
                levels = c("3", "4", "5", "6", "7", "8", "9", "10-14", "15-19", "\u2265 20")),
  Count = c(binCount_case, binCount_control),
  Group = rep(c("Cases (n = 1362)", "Controls (n = 1026)"),each = 10)
)

# Make plot
p <- ggplot(plotDF) +
  geom_bar(aes(x = Name, fill = Group, y = Count), color = "white",
           stat = "identity", position = position_dodge()) +
  #facet_grid(rows = vars(Group)) +
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2') +
  scale_fill_manual(values = colors) +
  labs(fill = NULL) +
  xlab("# Samples") +
  ylab("# Datasets") +
  theme_minimal() +
  theme(legend.position = c(0.8,0.8),
        strip.text = element_blank())

# Save plot
ggsave(p, file = "3. PlotMeta/sampleBarchart1.png", width = 4.5, height = 3)


################################################################################

# Donut chart of selected diseases

################################################################################

# select diseases
selDis <- c("FXS", "DMD", "DS", "RTT")
metaData$Disease1 <- metaData$Disease
metaData$Disease1[!(metaData$`Disease abbreviation` %in% selDis)] <- "Other"


# Create test data.
data <- data.frame(
  category= names(table(metaData$Disease1)),
  count= as.numeric(table(metaData$Disease1))
)

data$category <- factor(data$category,
                        levels =  c("Rett Syndrome", "Duchenne Muscular Dystrophy", "Fragile X Syndrome",
                                    "Down Syndrome","Other"))

data <- arrange(data, by = category)

# Compute percentages
data$fraction = data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n=-1))

# Set colors
colors <- setNames(c("#D95F02", "#7570B3", "#E7298A", "#E6AB02", "#D9D9D9"),
         c("Rett Syndrome", "Duchenne Muscular Dystrophy", "Fragile X Syndrome",
            "Down Syndrome","Other"))

# Make the plot
p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) +
  labs(fill = NULL) +
  scale_fill_manual(values = colors) +
  theme_void() 

# Save plot
ggsave(p, file = "3. PlotMeta/DiseasePiechart.png", width = 6, height = 5)

