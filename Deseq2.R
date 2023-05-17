library(DESeq2)
library(ggplot2)
library(readr)
library(pheatmap)
library(RColorBrewer)
library (vsn)
library(dplyr)
library(gplots)
library(AnnotationDbi)
library (EnsDb.Mmusculus.v79)

### CREATE A DICTIONARY OF GENEID FROM ENSEMBL DATABASE
keys <- keys(EnsDb.Mmusculus.v79, keytype="GENEID")
maps <- select(EnsDb.Mmusculus.v79, keys=keys, columns=c("GENENAME", "GENEBIOTYPE"),
               keytype="GENEID")

#All data
setwd()
countData <- read.csv(counts.csv", header = TRUE, row.names = 1)
countData <-as.matrix(countData)
head(countData)
metaData <- read.csv("meta.csv", header = TRUE)
metaData

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~genotype)
dds

#prefilter low count genes
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]
dds
#assign treatment levels.  This will make cisplatin the numerator and saline the denominator
dds$genotype <- factor(dds$genotype, levels = c("female_dehy", "female_ad", "male_dehy", "male_ad"))

#DEG
dds <-DESeq(dds)

result <-results(dds, contrast=c("genotype", "female_dehy", "female_ad"), alpha = 0.05)
head(result)
summary(result, )
normalized_counts <- counts(dds, normalized=TRUE)

plotMA(result, main=paste0('Condition: Control vs. Treatment'), ylim=c(-5,5))
rld <- rlogTransformation(dds, blind = TRUE)
plotPCA(rld, intgroup = c("genotype"))

#Plot counts for a single gene. Below is the plot for the gene with the lowest p-value:
plotCounts(dds, gene=which.min(result$padj), intgroup='genotype', pch = 19)

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
vsd <- vst(dds, blind=FALSE)

meanSdPlot(assay(vsd))

#heatmap of samples
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype", "ID")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap and sample to sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$id, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#look for sample outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

## Merge with normalized count data
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"

## Annotate Ensembl ID using Ensembl database EnsDb.Mmusculus.v79
resdata$GeneName <- with(resdata, maps$GENENAME[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)
resdata$GeneType <- with(resdata, maps$GENEBIOTYPE[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="midkiodf_results.csv")

#normalized counts
normalized_counts <- merge(as.data.frame(normalized_counts), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(normalized_counts)[1] <- "EnsemblID"

## Annotate Ensembl ID using Ensembl database EnsDb.Mmusculus.v79
normalized_counts$GeneName <- with(normalized_counts, maps$GENENAME[match(EnsemblID, maps$GENEID)])
normalized_counts <- normalized_counts %>% relocate(GeneName, .after = EnsemblID)
normalized_counts$GeneType <- with(normalized_counts, maps$GENEBIOTYPE[match(EnsemblID, maps$GENEID)])
normalized_counts <- normalized_counts %>% relocate(GeneType, .after = GeneName)

write.table(normalized_counts, file="normalized_countsf.txt", sep="\t", quote=F, col.names=NA)


#################Male dehydrated vs ad#############################################
result <-results(dds, contrast=c("genotype", "male_dehy", "male_ad"), alpha = 0.05)
head(result)
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"
## Annotate Ensembl ID using Ensembl database EnsDb.Mmusculus.v79
resdata$GeneName <- with(resdata, maps$GENENAME[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)
resdata$GeneType <- with(resdata, maps$GENEBIOTYPE[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="midkiodm_results.csv")


####################Male vs female ad##############################################
result <-results(dds, contrast=c("genotype", "male_ad", "female_ad"), alpha = 0.05)
head(result)
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"
## Annotate Ensembl ID using Ensembl database EnsDb.Mmusculus.v79
resdata$GeneName <- with(resdata, maps$GENENAME[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)
resdata$GeneType <- with(resdata, maps$GENEBIOTYPE[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="midkidMvFadlib_results.csv")

###################Male vs female dehydrated########################################
result <-results(dds, contrast=c("genotype", "male_dehy", "female_dehy"), alpha = 0.05)
head(result)
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"
## Annotate Ensembl ID using Ensembl database EnsDb.Mmusculus.v79
resdata$GeneName <- with(resdata, maps$GENENAME[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)
resdata$GeneType <- with(resdata, maps$GENEBIOTYPE[match(EnsemblID, maps$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="midkidMvFdehydrated_results.csv")

#####################Making Venn Diagrams###############################
#need to take the numbers calculated for unique and overlap and write the code in a way that it will plot it properly.  In this case we have 2219 common genes, 629 genes in ad lib and 2 genes in dehydrated for greater in males.  For colors go here https://venn.bio-spring.top/using-ggvenndiagram

library(ggVennDiagram)
library(ggplot2)

#male data
x <- list(
  Ad_lib = 3:2850, 
  Dehydrated = 1:2221
)
p <-ggVennDiagram(x)
# Blues for the males
p + scale_fill_distiller(palette = "Blues", direction = 1)

#Female data
x <- list(
  Ad_lib = 1:2596, 
  Dehydrated = 909:4619
)
p2 <-ggVennDiagram(x)
# Reds for the females
p2 + scale_fill_distiller(palette = "Reds", direction = 1)

p+p2+ scale_fill_distiller(palette = "Reds", direction = 1)

