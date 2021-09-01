##### RPK cutoff determination #####
library(reshape2)
library(ggplot2)
library(ggpubr)

setwd("/path/to/working-directory/")

ftable <- read.csv("rpkNonMaskedDirtyTable.txt", row.names=1, header=TRUE, sep='\t')

singleColumn <- as.vector(as.matrix(ftable[,colnames(ftable)]))
noZeros <- singleColumn[singleColumn!=0]

#hist(noZeros)
frequencyDistribution <- table(noZeros)
quantile(noZeros, .1)

richnessTable = ftable[FALSE, ]
for(i in seq(0, 1, by = 0.01)) {
  richnessTable <- rbind(richnessTable, colSums(ftable > quantile(noZeros, i)))
}
row.names(richnessTable) <- seq(0, 100, by=1)
colnames(richnessTable) <- colnames(ftable)

finalTable <- melt(as.matrix(richnessTable))
colnames(finalTable) <- c("Percentile", "Sample", "Richness")

p1 <- ggplot(data = finalTable, aes(x=Percentile, y=Richness)) + 
  geom_line(aes(colour=Sample)) +
  scale_x_reverse() + 
  theme(aspect.ratio = 1) +
  ggtitle("100 to 0 percentile")

# Zoom in
richnessTable = ftable[FALSE, ]
for(i in seq(0, 0.25, by = 0.01)) {
  richnessTable <- rbind(richnessTable, colSums(ftable > quantile(noZeros, i)))
}
row.names(richnessTable) <- seq(0, 25, by=1)
colnames(richnessTable) <- colnames(ftable)

finalTable <- melt(as.matrix(richnessTable))
colnames(finalTable) <- c("Percentile", "Sample", "Richness")

p2 <- ggplot(data = finalTable, aes(x=Percentile, y=Richness)) + 
  geom_line(aes(colour=Sample)) +
  scale_x_reverse() + 
  theme(aspect.ratio = 1) +
  ggtitle("25 to 0 percentile")

# Zoom in in
richnessTable = ftable[FALSE, ]
for(i in seq(0, 0.1, by = 0.01)) {
  richnessTable <- rbind(richnessTable, colSums(ftable > quantile(noZeros, i)))
}
row.names(richnessTable) <- seq(0, 10, by=1)
colnames(richnessTable) <- colnames(ftable)

finalTable <- melt(as.matrix(richnessTable))
colnames(finalTable) <- c("Percentile", "Sample", "Richness")

p3 <- ggplot(data = finalTable, aes(x=Percentile, y=Richness)) + 
  geom_line(aes(colour=Sample)) +
  scale_x_reverse() + 
  theme(aspect.ratio = 1) +
  ggtitle("10 to 0 percentile")

ggarrange(p1, p2, p3, nrow = 1, 
          common.legend=TRUE, 
          legend = "bottom")

##### Mask table #####
setwd("/path/to/working-directory/")

ftable <- read.csv("rpkNonMaskedDirtyTable.txt", row.names=1, header=TRUE, sep='\t')

ftable[ftable < 1] <- 0
maskedTable <- ftable[rowSums(ftable) > 0,]

write.table(maskedTable, "rpkMaskedDirty.txt", sep='\t', col.names=NA)

##### Decontam ######
library(decontam)
library(phyloseq)
library(gridExtra)

setwd("/path/to/working-directory/")

ftable <- read.csv("rpkMaskedDirty.txt", row.names=1, header=TRUE, sep='\t')
data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]
data$Lactobacillus.dominance...80.. <- as.factor(data$Lactobacillus.dominance...80..)

frequency_table <- as.matrix(t(ftable))
contam <- isContaminant(frequency_table, method='prevalence', neg=data$Batch =='Buffer', threshold=0.1)
clean_frequency_table <- ftable[!contam$contaminant=="TRUE", ]
clean_frequency_table <- clean_frequency_table[rowSums(clean_frequency_table) > 0,]

presenceAbsence <- clean_frequency_table
presenceAbsence[presenceAbsence > 0] <- 1
presenceAbsence <- presenceAbsence[,colSums(presenceAbsence) > 0]

OTU = otu_table(presenceAbsence, taxa_are_rows = TRUE)
metadata = sample_data(data)
physeq = phyloseq(OTU, metadata)
dist.ord <- ordinate(physeq, "PCoA", "bray")
lacDom = plot_ordination(physeq, dist.ord, color="Lactobacillus.dominance...80..", title="LacDom") + theme(aspect.ratio =1)
giScore = plot_ordination(physeq, dist.ord, color="Binary.genital.inflammatory.score", title="GI Score") + theme(aspect.ratio =1)
grid.arrange(lacDom, giScore, ncol=1)

table(contam$contaminant)
colSums(clean_frequency_table)
buffersRemoved <- subset(clean_frequency_table, select=-c(I2058, I2073))

write.table(buffersRemoved, "rpkMaskedClean.txt", sep='\t', col.names=NA)

##### Diversity metrics ######
library(phyloseq)
library(vegan)
library(ggpubr)
library(gridExtra)
library(ggplot2)
library(rich)

setwd("/path/to/working-directory/")
ftable <- read.csv("rpkMaskedClean.txt", row.names=1, header=TRUE, sep='\t')

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]
data$Lactobacillus.dominance...80.. <- as.factor(data$Lactobacillus.dominance...80..)

# Bray Curtis PCoAs
presenceAbsence <- ftable
presenceAbsence[presenceAbsence > 0] <- 1
OTU = otu_table(presenceAbsence, taxa_are_rows = TRUE)
metadata = sample_data(data)
physeqPresenceAbsence = phyloseq(OTU, metadata)
dist.ord <- ordinate(physeqPresenceAbsence, "PCoA", "bray")
lacDom = plot_ordination(physeqPresenceAbsence, dist.ord, color="Lactobacillus.dominance...80..", title="LacDom") + theme(aspect.ratio =1)
giScore = plot_ordination(physeqPresenceAbsence, dist.ord, color="Binary.genital.inflammatory.score", title="GI Score") + theme(aspect.ratio =1)
grid.arrange(lacDom, giScore, ncol=1, top="All viruses contig RPK")

adonisLac <- as.numeric(data$Lactobacillus.dominance...80..)
adon.results<-adonis(t(ftable) ~ adonisLac, method="bray",perm=10000)
print(adon.results)

adonisGI <- data$Binary.genital.inflammatory.score
adon.results<-adonis(t(ftable) ~ adonisGI, method="bray",perm=10000)
print(adon.results)

results <- data.frame(row.names=row.names(data),
                      Lactobacillus=data$Lactobacillus.dominance...80..,
                      GIScore=data$Binary.genital.inflammatory.score,
                      PC1=lacDom$data$Axis.1,
                      PC2=lacDom$data$Axis.2)

# Richness/Shannon
results$Shannon <- diversity(t(ftable), index="shannon")
results$Richness <- colSums(ftable != 0)

write.table(results, "diversityResults.txt", sep='\t', col.names=NA)

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
grid.arrange(p1,p2,p3,p4, nrow=2, top="All viruses contig RPK")

# Richness rarefaction curve
#noInflam <- ftable[,data$Binary.genital.inflammatory.score=="0"]
#inflam <- ftable[,data$Binary.genital.inflammatory.score==">0"]
#rarcNoInflam <- rarc(t(noInflam), samplesize=c(0:15), nrandom=100)
#rarcInflam <- rarc(t(inflam), samplesize=c(0:8), nrandom=100)

# Relative abundance
tax <- read.csv("contigTaxonomy.txt", row.names = 1, header = TRUE, sep = '\t')
tax <- as.matrix(tax[row.names(ftable),])

physeq = phyloseq(otu_table(ftable, taxa_are_rows = TRUE),
                  sample_data(data),
                  tax_table(tax))

relAbundance = transform_sample_counts(physeq, function(x) {x / sum(x)})
relPlot <- plot_bar(relAbundance, fill="family")
relPlot + geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

collapsedFamily <- tax_glom(relAbundance, "family")
relAbundanceExport <- as(otu_table(collapsedFamily), "matrix")
familyRelativeAbundance <- as.data.frame(relAbundanceExport)
familyTaxonomy <- tax[row.names(familyRelativeAbundance),"family"]
row.names(familyRelativeAbundance) <- familyTaxonomy

write.table(t(familyRelativeAbundance), "familyRelativeAbundanceCounts.txt", sep='\t', col.names=NA)

##### LefSe heatmap #####
library(gplots)
library(RColorBrewer)

setwd("/path/to/working-directory/")
ftable <- read.csv("rpkMaskedClean.txt", row.names=1, header=TRUE, sep='\t')

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]

tax <- read.csv("contigTaxonomy.txt", row.names = 1, header=TRUE, sep = '\t')

header <- c("Contig", "tempLDA", "meta", "LDA", "p")
lacLDA <- read.csv("lefse/lacLDA.txt", header=FALSE, sep='\t')
colnames(lacLDA) <- header
lacLDA <- lacLDA[!is.na(lacLDA$LDA),]

giLDA <- read.csv("lefse/giLDA.txt", header=FALSE, sep='\t')
colnames(giLDA) <- header
giLDA <- giLDA[!is.na(giLDA$LDA),]

lacFtable <- ftable[lacLDA$Contig, ]
giFtable <- ftable[giLDA$Contig, ]

leftLacTable <- lacFtable[, data$Lactobacillus.dominance...80.. == "1"]
rightLacTable <- lacFtable[, data$Lactobacillus.dominance...80.. == "0"]

leftGITable <- giFtable[, data$Binary.genital.inflammatory.score == "0"]
rightGITable <- giFtable[, data$Binary.genital.inflammatory.score == ">0"]

leftLacOrder <- hclust(dist(t(leftLacTable)))$order
rightLacOrder <- hclust(dist(t(rightLacTable)))$order

leftGIOrder <- hclust(dist(t(leftGITable)))$order
rightGIOrder <- hclust(dist(t(rightGITable)))$order

leftLacTable <- leftLacTable[, leftLacOrder]
rightLacTable <- rightLacTable[, rightLacOrder]

leftGITable <- leftGITable[, leftGIOrder]
rightGITable <- rightGITable[, rightGIOrder]

orderedLacTable <- cbind(leftLacTable, rightLacTable)
orderedGITable <- cbind(leftGITable, rightGITable)

orderedLacTable$family <- tax[row.names(orderedLacTable), "family"]
orderedGITable$family <- tax[row.names(orderedGITable), "family"]

topLacTable <- orderedLacTable[lacLDA$meta == 0, ]
bottomLacTable <- orderedLacTable[lacLDA$meta == 1, ]

topGITable <- orderedGITable[giLDA$meta == ">0", ]
bottomGITable <- orderedGITable[giLDA$meta == "0", ]

topLacTable$family <- factor(topLacTable$family, 
                             levels = names(sort(table(topLacTable$family), decreasing=TRUE)))
topLacTable <- topLacTable[order(topLacTable$family), ]

bottomLacTable$family <- factor(bottomLacTable$family, 
                             levels = names(sort(table(bottomLacTable$family), decreasing=TRUE)))
bottomLacTable <- bottomLacTable[order(bottomLacTable$family), ]

topGITable$family <- factor(topGITable$family, 
                                levels = names(sort(table(topGITable$family), decreasing=TRUE)))
topGITable <- topGITable[order(topGITable$family), ]

bottomGITable$family <- factor(bottomGITable$family, 
                                levels = names(sort(table(bottomGITable$family), decreasing=TRUE)))
bottomGITable <- bottomGITable[order(bottomGITable$family), ]

finalLacTable <- rbind(topLacTable, bottomLacTable)
finalLacTable$family <- NULL

finalGITable <- rbind(topGITable, bottomGITable)
finalGITable$family <- NULL

logLacTable <- log10(data.matrix(finalLacTable)+1)
logGITable <- log10(data.matrix(finalGITable)+1)

# Add family colors in
lacFamily <- factor(tax[row.names(logLacTable), "family"])
lacNumFamilies <- length(unique(lacFamily))
lacRowColors <- brewer.pal(lacNumFamilies, "Set1")

giFamily <- factor(tax[row.names(logGITable), "family"])
giNumFamilies <- length(unique(lacFamily))
giRowColors <- brewer.pal(giNumFamilies, "Set1")

#giOrder <- data[order(data$Genital.inflammatory.score),]
giOrder <- data[data$Genital.inflammatory.score,]
giScore <- factor(giOrder$Genital.inflammatory.score)
logGITable <- logGITable[, row.names(giOrder)]
giScoreColors <- brewer.pal(length(unique(giScore)), "Oranges")
giScoreColors[1] <- "#FFFFFF"

heatmap.2(logLacTable, col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=14, rowsep=13, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Lactobacillus LEfSe heatmap Log10",
          RowSideColors=lacRowColors[lacFamily])

legend(x="left", legend=levels(lacFamily),
       fill=lacRowColors)

heatmap.2(logGITable, col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=15, rowsep=40, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Inflammation LEfSe heatmap",
          RowSideColors=giRowColors[giFamily],
          ColSideColors=giScoreColors[giScore])

# Unsupervised clustering
heatmap.2(logGITable, col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=15, rowsep=40, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Inflammation LEfSe heatmap",
          RowSideColors=giRowColors[giFamily])

legend(x="left", legend=levels(giFamily),
       fill=giRowColors)
