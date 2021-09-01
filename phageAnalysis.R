##### RPK cutoff determination #####
library(reshape2)
library(ggplot2)
library(ggpubr)

setwd("/path/to/working-directory/")

ftable <- read.csv("phages/phageRPKNonMaskedDirtyTable.txt", row.names=1, header=TRUE, sep='\t')

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

ftable <- read.csv("phages/phageRPKNonMaskedDirtyTable.txt", row.names=1, header=TRUE, sep='\t')

ftable[ftable < 1] <- 0
maskedTable <- ftable[rowSums(ftable) > 0,]

write.table(maskedTable, "phages/phageRpkMaskedDirty.txt", sep='\t', col.names=NA)

##### Decontam ######
library(decontam)
library(phyloseq)
library(gridExtra)
library(ggplot2)

setwd("/path/to/working-directory/")

ftable <- read.csv("phages/phageRpkMaskedDirty.txt", row.names=1, header=TRUE, sep='\t')
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

write.table(clean_frequency_table, "phages/phageRpkMaskedCleanWithBuffer.txt", sep='\t', col.names=NA)
write.table(buffersRemoved, "phages/phageRpkMaskedClean.txt", sep='\t', col.names=NA)

##### Diversity metrics ######
library(phyloseq)
library(vegan)
library(ggpubr)
library(gridExtra)
library(ggplot2)
library(rich)

setwd("/path/to/working-directory/")
ftable <- read.csv("phages/phageRpkMaskedClean.txt", row.names=1, header=TRUE, sep='\t')

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
grid.arrange(lacDom, giScore, ncol=1, top="Bacteriophages contig RPK")

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

write.table(results, "phages/phageDiversityResults.txt", sep='\t', col.names=NA)

p1 <- ggboxplot(results, x="Lactobacillus", y="Shannon", 
                title="Phage Lactobacillus Shannon", xlab="Lactobacillus dominance") +
  stat_compare_means(comparisons=list(c("0", "1")))
p2 <- ggboxplot(results, x="Lactobacillus", y="Richness", 
                title="Phage Lactobacillus Richness", xlab="Lactobacillus dominance") +
  stat_compare_means(comparisons=list(c("0", "1")))
p3 <- ggboxplot(results, x="GIScore", y="Shannon", 
                title="Phage genital inflammation Shannon", xlab="Genital inflammation score") + 
  stat_compare_means(comparisons=list(c("0", ">0")))
p4 <- ggboxplot(results, x="GIScore", y="Richness",
                title="Phage genital inflammation Richness", xlab="Genital inflammation score") +
  stat_compare_means(comparisons=list(c("0", ">0")))
grid.arrange(p1,p2,p3,p4, nrow=2, top="Bacteriophages contig RPK")


noInflam <- presenceAbsence[,data$Binary.genital.inflammatory.score=="0"]
inflam <- presenceAbsence[,data$Binary.genital.inflammatory.score==">0"]
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

write.table(t(familyRelativeAbundance), "phages/phageFamilyRelativeAbundanceCounts.txt", sep='\t', col.names=NA)

##### LefSe heatmap #####
library(gplots)
library(RColorBrewer)

setwd("/path/to/working-directory/")
ftable <- read.csv("phages/phageRpkMaskedClean.txt", row.names=1, header=TRUE, sep='\t')

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]

tax <- read.csv("contigTaxonomy.txt", row.names = 1, header=TRUE, sep = '\t')

header <- c("Contig", "tempLDA", "meta", "LDA", "p")
lacLDA <- read.csv("phages/lefse/lacLDA.txt", header=FALSE, sep='\t')
colnames(lacLDA) <- header
lacLDA <- lacLDA[!is.na(lacLDA$LDA),]

giLDA <- read.csv("phages/lefse/giLDA.txt", header=FALSE, sep='\t')
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

heatmap.2(logLacTable, col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=14, rowsep=9, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Phage Lactobacillus LEfSe heatmap (Log)",
          RowSideColors=lacRowColors[lacFamily])

legend(x="left", legend=levels(lacFamily),
       fill=lacRowColors)

heatmap.2(logGITable, col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=15, rowsep=18, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Inflammation LEfSe heatmap (Log)",
          RowSideColors=giRowColors[giFamily])

legend(x="left", legend=levels(giFamily),
       fill=giRowColors)

# No log transform
lacFamily <- factor(tax[row.names(finalLacTable), "family"])
lacNumFamilies <- length(unique(lacFamily))
lacRowColors <- brewer.pal(lacNumFamilies, "Set1")

giFamily <- factor(tax[row.names(finalGITable), "family"])
giNumFamilies <- length(unique(lacFamily))
giRowColors <- brewer.pal(giNumFamilies, "Set1")

heatmap.2(data.matrix(finalLacTable), col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=14, rowsep=219, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Lactobacillus LEfSe heatmap",
          RowSideColors=lacRowColors[lacFamily])

legend(x="left", legend=levels(lacFamily),
       fill=lacRowColors)

heatmap.2(data.matrix(finalGITable), col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE,
          colsep=15, rowsep=18, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Inflammation LEfSe heatmap",
          RowSideColors=giRowColors[giFamily])

legend(x="left", legend=levels(giFamily),
       fill=giRowColors)

write.table(finalLacTable, "phages/lefse/lacHeatmapTable.txt", sep='\t', col.names=NA)
write.table(finalGITable, "phages/lefse/giHeatmapTable.txt", sep='\t', col.names=NA)

# Strictly hierarchical clustering
library(gplots)
library(RColorBrewer)

setwd("/path/to/working-directory/")
ftable <- read.csv("phages/phageRpkMaskedClean.txt", row.names=1, header=TRUE, sep='\t')

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]

tax <- read.csv("contigTaxonomy.txt", row.names = 1, header=TRUE, sep = '\t')

header <- c("Contig", "tempLDA", "meta", "LDA", "p")
lacLDA <- read.csv("phages/lefse/lacLDA.txt", header=FALSE, sep='\t')
colnames(lacLDA) <- header
lacLDA <- lacLDA[!is.na(lacLDA$LDA),]

lacFtable <- ftable[lacLDA$Contig, ]

leftLacTable <- lacFtable[, data$Lactobacillus.dominance...80.. == "1"]
rightLacTable <- lacFtable[, data$Lactobacillus.dominance...80.. == "0"]

leftLacOrder <- hclust(dist(t(leftLacTable)))$order
rightLacOrder <- hclust(dist(t(rightLacTable)))$order

leftLacTable <- leftLacTable[, leftLacOrder]
rightLacTable <- rightLacTable[, rightLacOrder]

orderedLacTable <- cbind(leftLacTable, rightLacTable)

orderedLacTable$family <- tax[row.names(orderedLacTable), "family"]

topLacTable <- orderedLacTable[lacLDA$meta == 0, ]
bottomLacTable <- orderedLacTable[lacLDA$meta == 1, ]

topLacTable$family <- factor(topLacTable$family, 
                             levels = names(sort(table(topLacTable$family), decreasing=TRUE)))
topLacTable <- topLacTable[order(topLacTable$family), ]

bottomLacTable$family <- factor(bottomLacTable$family, 
                                levels = names(sort(table(bottomLacTable$family), decreasing=TRUE)))
bottomLacTable <- bottomLacTable[order(bottomLacTable$family), ]

finalLacTable <- rbind(topLacTable, bottomLacTable)
finalLacTable$family <- NULL

logLacTable <- log10(data.matrix(finalLacTable)+1)

lacPercentages <- data[colnames(logLacTable), "Lactobacillus.abundance"]
colors <- colorRampPalette(c('white', 'blue'))
colorRange <- colors(101)
lacColors <- c()
for(val in lacPercentages)
  lacColors <- c(lacColors, colorRange[round(val*100)+1])

# Add family colors in
lacFamily <- factor(tax[row.names(logLacTable), "family"])
lacNumFamilies <- length(unique(lacFamily))
lacRowColors <- brewer.pal(lacNumFamilies, "Set1")

heatmap.2(logLacTable, col=brewer.pal(9,"Blues"),
          Rowv=FALSE, Colv=FALSE, 
          rowsep=9, sepcolor="Black",
          trace="none", density.info="none",
          labCol=FALSE, labRow=FALSE, main="Phage Lactobacillus LEfSe heatmap (Log) Clustering",
          RowSideColors=lacRowColors[lacFamily],
          ColSideColors=lacColors)

legend(x="left", legend=levels(lacFamily),
       fill=lacRowColors)

##### Transkingdom analysis #####
library(corrplot)
library(gplots)

setwd("/path/to/working-directory/")
ftable <- read.csv("phages//bacteriaPhageInteraction/correlation.txt", row.names=1, header=TRUE, sep='\t')
ftable <- ftable[rowSums(ftable) != 0,]

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]

noInflam <- t(ftable[, row.names(data[data$Binary.genital.inflammatory.score==0,])])
noInflam <- noInflam[, colSums(noInflam) != 0]

noInflamRows <- colnames(noInflam)

noInflamCor <- cor(noInflam, method="pearson")

noInflamCor <- noInflamCor[-grep("Contig", row.names(noInflamCor)),
                           grep("Contig", row.names(noInflamCor))]

noInflamHeatmap <- heatmap.2(noInflamCor,
                             col = redblue(100),
                             trace = "none",
                             density.info = "none",
                             xlab = "Bacteriophages",
                             ylab = "Bacteria",
                             main = "No inflammation",
                             labRow = "",
                             labCol = "")

noInflamHeatmap

inflam <- t(ftable[noInflamRows, row.names(data[data$Binary.genital.inflammatory.score==">0",])])
inflamCor <- cor(inflam, method="pearson")

inflamCor <- inflamCor[-grep("Contig", row.names(inflamCor)),
                           grep("Contig", row.names(inflamCor))]

inflamCor <- inflamCor[noInflamHeatmap$rowInd, noInflamHeatmap$colInd]

inflamHeatmap <- heatmap.2(inflamCor,
                           col = redblue(100),
                           trace = "none",
                           density.info = "none",
                           Rowv = FALSE,
                           Colv = FALSE,
                           dendrogram = "none",
                           xlab = "Bacteriophages",
                           ylab = "Bacteria",
                           main = "Inflammation",
                           labRow = "",
                           labCol = "")
inflamHeatmap

sortedNoInflamCor <- noInflamCor[rev(noInflamHeatmap$rowInd), noInflamHeatmap$colInd]
sortedInflamCor <- inflamCor[rev(row.names(inflamCor)),]

write.table(sortedNoInflamCor, "phages/bacteriaPhageInteraction/sortedNoInflamCor.txt", sep='\t', col.names=NA)
write.table(inflamCor, "phages/bacteriaPhageInteraction/sortedInflamCor.txt", sep='\t', col.names=NA)

hc <- as.hclust(noInflamHeatmap$colDendrogram)
right <- cutree(hc, k=4)
table(right)

hc <- as.hclust(noInflamHeatmap$rowDendrogram)
right <- cutree(hc, k=4)
table(right)

### Transkingom clustered by inflammation
library(corrplot)
library(gplots)
library(ggplot2)
library(reshape2)
library(ggforce)

setwd("/path/to/working-directory/")
ftable <- read.csv("phages//bacteriaPhageInteraction/correlation.txt", row.names=1, header=TRUE, sep='\t')
ftable <- ftable[rowSums(ftable) != 0,]

data <- read.csv("CVLMetadata.txt", row.names = 1, header=TRUE, sep = '\t')
data <- data[colnames(ftable),]

inflam <- t(ftable[, row.names(data[data$Binary.genital.inflammatory.score==">0",])])
inflam <- inflam[, colSums(inflam) != 0]

inflamRows <- colnames(inflam)

inflamCor <- cor(inflam, method="pearson")

inflamCor <- inflamCor[-grep("Contig", row.names(inflamCor)),
                       grep("Contig", row.names(inflamCor))]

inflamHeatmap <- heatmap.2(inflamCor,
                           col = redblue(100),
                           trace = "none",
                           density.info = "none",
                           xlab = "Bacteriophages",
                           ylab = "Bacteria",
                           main = "Inflammation",
                           labRow = "",
                           labCol = "")
inflamHeatmap

# No inflammation in same row/column order as inflammation
noInflam <- t(ftable[inflamRows, row.names(data[data$Binary.genital.inflammatory.score==0,])])
noInflamCor <- cor(noInflam, method="pearson")

noInflamCor <- noInflamCor[-grep("Contig", row.names(noInflamCor)),
                           grep("Contig", row.names(noInflamCor))]

noInflamCor <- noInflamCor[inflamHeatmap$rowInd, inflamHeatmap$colInd]

noInflamHeatmap <- heatmap.2(noInflamCor,
                             col = redblue(100),
                             trace = "none",
                             density.info = "none",
                             Rowv = FALSE,
                             Colv = FALSE,
                             dendrogram = "none",
                             xlab = "Bacteriophages",
                             ylab = "Bacteria",
                             main = "No inflammation",
                             labRow = "",
                             labCol = "")

noInflamHeatmap

sortedInflamCor <- inflamCor[rev(inflamHeatmap$rowInd), inflamHeatmap$colInd]
sortedNoInflamCor <- noInflamCor[rev(row.names(noInflamHeatmap)), noInflamHeatmap$colInd]

write.table(sortedInflamCor, "phages/bacteriaPhageInteraction/inflammationBased/sortedInflamCor.txt", sep='\t', col.names=NA)
write.table(sortedNoInflamCor, "phages/bacteriaPhageInteraction/inflammationBased/sortedNoInflamCor.txt", sep='\t', col.names=NA)

hc <- as.hclust(inflamHeatmap$colDendrogram)
table(cutree(hc, k=4))

hc <- as.hclust(inflamHeatmap$rowDendrogram)
table(cutree(hc, k=5))

# Dot plot of correlations between pathogenic bacteria and Lactobacillus
# Lactobacillus x-axis, pathogenic bacteria y-axis

# Subset rows from inflammation heatmap for 
# Mycoplasma, Prevotella, Sneathia, Mobiluncus, Gardnerella, Lactobacillus
mycoplasma <- "t__14984"
prevotella <- "t__12427"
sneathia <- "dnOTU_3"
mobiluncus <- "dnOTU_27"
gardnerella <- "t__91260"
lactobacillus <- "dnOTU_1"

# Pull rows out as vectors
lacto <- inflamCor[lactobacillus,]

df1 <- data.frame(x=lacto, y=inflamCor[mycoplasma,])
df2 <- data.frame(x=lacto, y=inflamCor[prevotella,])
df3 <- data.frame(x=lacto, y=inflamCor[sneathia,])
df4 <- data.frame(x=lacto, y=inflamCor[mobiluncus,])
df5 <- data.frame(x=lacto, y=inflamCor[gardnerella,])
final <- rbind(df1, df2, df3, df4, df5)

# Plot scatter plot
ggplot(final, aes(x,y)) + 
  geom_point(position=position_jitternormal(sd_x=0.025, sd_y=0.025), alpha=0.15, size=0.7) + 
  xlim(-1.1,1.1) + ylim(-1.1,1.1) + 
  xlab("Lactobacillus correlation") + ylab("Disease progression associated bacteria correlation") +
  ggtitle("Bacteriophage corrleation between \nLactobacillus and disease progression associated bacteria") +
  coord_fixed() +
  theme_light() +
  geom_smooth(method = lm, se=FALSE)

summary(lm(y~x, final))

# Anaerococcus chosen as control bacteria
bacteroides <- "dnOTU_57"
bacter <- inflamCor[bacteroides,]

df6 <- data.frame(x=bacter, y=inflamCor[mycoplasma,])
df7 <- data.frame(x=bacter, y=inflamCor[prevotella,])
df8 <- data.frame(x=bacter, y=inflamCor[sneathia,])
df9 <- data.frame(x=bacter, y=inflamCor[mobiluncus,])
df10 <- data.frame(x=bacter, y=inflamCor[gardnerella,])
controlFinal <- rbind(df6, df7, df8, df9, df10) 

ggplot(controlFinal, aes(x,y)) + 
  geom_point(position=position_jitternormal(sd_x=0.025, sd_y=0.025), alpha=0.15, size=0.7) + 
  xlim(-1.1,1.1) + ylim(-1.1,1.1) + 
  xlab("Bacteroides correlation") + ylab("Disease progression associated bacteria correlation") +
  ggtitle("Bacteriophage corrleation between \nBacteroides and disease progression associated bacteria") +
  coord_fixed() +
  theme_light() +
  geom_smooth(method = lm, se = FALSE)

summary(lm(y~x, controlFinal))
