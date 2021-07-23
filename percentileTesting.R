##### RPK cutoff determination #####
library(reshape2)
library(ggplot2)
library(ggpubr)
library(decontam)
library(phyloseq)
library(gridExtra)
library(vegan)
library(gridExtra)

START <- 0.1
END <- 0.9
STEP <- 0.1

setwd("E:\\Dropbox (ASU)\\Lim Lab\\Projects\\Cervicovaginal lavage (U of A)\\CVL\\Contigs\\allContigsAnalysis")

#ftable <- read.csv("rpkNonMaskedDirtyTable.txt", row.names=1, header=TRUE, sep='\t')
ftable <- read.csv("phages/phageRpkNonMaskedDirtyTable.txt", row.names=1, header=TRUE, sep='\t')

singleColumn <- as.vector(as.matrix(ftable[,colnames(ftable)]))
noZeros <- singleColumn[singleColumn!=0]

frequencyDistribution <- table(noZeros)

richnessTable = ftable[FALSE, ]
percentileCutoffs <- c()

for(i in seq(START, END, by = STEP)) {
  percentileCutoffs <- c(percentileCutoffs, quantile(noZeros, i))
  richnessTable <- rbind(richnessTable, colSums(ftable > quantile(noZeros, i)))
}
row.names(richnessTable) <- seq(START, END, by=STEP)
colnames(richnessTable) <- colnames(ftable)

finalTable <- melt(as.matrix(richnessTable))
colnames(finalTable) <- c("Percentile", "Sample", "Richness")

p1 <- ggplot(data = finalTable, aes(x=Percentile, y=Richness)) + 
  geom_line(aes(colour=Sample)) +
  scale_x_reverse() + 
  theme(aspect.ratio = 1) +
  ggtitle("10 to 25 percentile")

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')

percentile <- START
gridVec <- c()

for(cutoff in percentileCutoffs) {
  
  # Mask table
  maskedTable <- ftable
  maskedTable[maskedTable < cutoff] <- 0
  maskedTable <- maskedTable[rowSums(maskedTable) > 0,]
  
  # Decontam table
  frequency_table <- as.matrix(t(maskedTable))
  contam <- isContaminant(frequency_table, method='prevalence', neg=data$Batch =='Buffer', threshold=0.1)
  clean_frequency_table <- maskedTable[!contam$contaminant=="TRUE", ]
  clean_frequency_table <- clean_frequency_table[rowSums(clean_frequency_table) > 0,]
  buffersRemoved <- subset(clean_frequency_table, select=-c(I2058, I2073))
  
  # Diversity
  data <- data[colnames(buffersRemoved),]
  results <- data.frame(row.names=row.names(data),
                        Lactobacillus=data$Lactobacillus.dominance...80..,
                        GIScore=data$Binary.genital.inflammatory.score)
  
  results$Shannon <- diversity(t(buffersRemoved), index="shannon")
  results$Richness <- colSums(buffersRemoved != 0)
  
  p1 <- ggboxplot(results, x="Lactobacillus", y="Shannon", 
                  title="Lactobacillus Shannon", xlab="Lactobacillus dominance") +
    stat_compare_means(comparisons=list(c("0", "1")))
  p2 <- ggboxplot(results, x="Lactobacillus", y="Richness", 
                  title="Lactobacillus Richness", xlab="Lactobacillus dominance") +
    stat_compare_means(comparisons=list(c("0", "1")))
  p3 <- ggboxplot(results, x="GIScore", y="Shannon", 
                  title="Genital inflammation Shannon", xlab="Genital inflammation score") + 
    stat_compare_means(comparisons=list(c("0", ">0")))
  p4 <- ggboxplot(results, x="GIScore", y="Richness",
                  title="Genital inflammation Richness", xlab="Genital inflammation score") +
    stat_compare_means(comparisons=list(c("0", ">0")))
  
  gridVec <- c(gridVec, grid.arrange(p1,p2,p3,p4, nrow=2, top=paste("All viruses at", percentile)))
  percentile <- percentile + STEP
}
percentileCutoffs
