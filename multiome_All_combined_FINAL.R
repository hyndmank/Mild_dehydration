library(Seurat)#v4.3.0
library(ggplot2)
library(Signac)#v1.9
library(dplyr)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(patchwork)
set.seed(1234)                
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(GenomicRanges)
library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(limma)
library(openxlsx)
library(tidyverse)
(future.globals.maxSize = 150000 * 1024^2)


# create a Seurat object containing the RNA adata
# load the RNA data
counts.630 <- Read10X("~/Desktop/Dehydration/630/outs/filtered_feature_bc_matrix")
counts.636 <- Read10X("~/Desktop/Dehydration/636/outs/filtered_feature_bc_matrix")
counts.638 <- Read10X("~/Desktop/Dehydration/638/outs/filtered_feature_bc_matrix")
counts.641 <- Read10X("~/Desktop/Dehydration/641/outs/filtered_feature_bc_matrix")

# load metadata
metadata.630 <- read.table(
  file = "~/Desktop/Dehydration/630/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

metadata.636 <- read.table(
  file = "~/Desktop/Dehydration/636/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

metadata.638 <- read.table(
  file = "~/Desktop/Dehydration/638/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

metadata.641 <- read.table(
  file = "~/Desktop/Dehydration/641/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

#load the ATAC fragments
frags.630 <- CreateFragmentObject(
  path = "~/Desktop/Dehydration/630/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata.630))

frags.636 <- CreateFragmentObject(
  path = "~/Desktop/Dehydration/636/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata.636))

frags.638 <- CreateFragmentObject(
  path = "~/Desktop/Dehydration/638/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata.638))

frags.641 <- CreateFragmentObject(
  path = "~/Desktop/Dehydration/641/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata.641))

# get gene annotations for mm10
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# create a Seurat object containing the RNA adata
sample_630 <- CreateSeuratObject(
  counts = counts.630$`Gene Expression`,
  assay = "RNA",
  meta.data = metadata.630
)

sample_636 <- CreateSeuratObject(
  counts = counts.636$`Gene Expression`,
  assay = "RNA",
  meta.data = metadata.636
)

sample_638 <- CreateSeuratObject(
  counts = counts.638$`Gene Expression`,
  assay = "RNA",
  meta.data = metadata.638
)

sample_641 <- CreateSeuratObject(
  counts = counts.641$`Gene Expression`,
  assay = "RNA",
  meta.data = metadata.641
)

# create ATAC assay and add it to the object,the assya ATAC is the original peaks called by cellranger.
sample_630[["ATAC"]] <- CreateChromatinAssay(
  counts = counts.630$Peaks,
  sep = c(":", "-"),
  fragments = frags.630,
  min.cells = 10,
  annotation = annotations
)
sample_630
#An object of class Seurat 
#149395 features across 5993 samples within 2 assays 
#Active assay: RNA (32285 features, 0 variable features)
#1 other assay present: ATAC

sample_636[["ATAC"]] <- CreateChromatinAssay(
  counts = counts.636$Peaks,
  sep = c(":", "-"),
  fragments = frags.636,
  min.cells = 10,
  annotation = annotations
)
sample_636
#An object of class Seurat 
#146215 features across 4915 samples within 2 assays 
#Active assay: RNA (32285 features, 0 variable features)
#1 other assay present: ATAC

sample_638[["ATAC"]] <- CreateChromatinAssay(
  counts = counts.638$Peaks,
  sep = c(":", "-"),
  fragments = frags.638,
  min.cells = 10,
  annotation = annotations
)
sample_638
#An object of class Seurat 
#152665 features across 5775 samples within 2 assays 
#Active assay: RNA (32285 features, 0 variable features)
#1 other assay present: ATAC

sample_641[["ATAC"]] <- CreateChromatinAssay(
  counts = counts.641$Peaks,
  sep = c(":", "-"),
  fragments = frags.641,
  min.cells = 10,
  annotation = annotations
)
sample_641
#An object of class Seurat 
#149036 features across 6138 samples within 2 assays 
#Active assay: RNA (32285 features, 0 variable features)
#1 other assay present: ATAC

#Q&C
DefaultAssay(sample_630) <- "ATAC"
sample_630 <- NucleosomeSignal(sample_630)
sample_630$nucleosome_group <- ifelse(sample_630$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_630, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
sample_630 <- TSSEnrichment(sample_630, fast=FALSE)
sample_630$high.tss <- ifelse(sample_630$TSS.enrichment > 2, 'High', 'Low')

VlnPlot(
  object = sample_630,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

DefaultAssay(sample_636) <- "ATAC"
sample_636 <- NucleosomeSignal(sample_636)
sample_636$nucleosome_group <- ifelse(sample_636$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_636, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
sample_636 <- TSSEnrichment(sample_636, fast=FALSE)
sample_636$high.tss <- ifelse(sample_636$TSS.enrichment > 2, 'High', 'Low')

VlnPlot(
  object = sample_636,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

DefaultAssay(sample_638) <- "ATAC"
sample_638 <- NucleosomeSignal(sample_638)
sample_638$nucleosome_group <- ifelse(sample_638$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_638, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
sample_638 <- TSSEnrichment(sample_638, fast=FALSE)
sample_638$high.tss <- ifelse(sample_638$TSS.enrichment > 2, 'High', 'Low')

VlnPlot(
  object = sample_638,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

DefaultAssay(sample_641) <- "ATAC"
sample_641 <- NucleosomeSignal(sample_641)
sample_641$nucleosome_group <- ifelse(sample_641$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_641, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
sample_641 <- TSSEnrichment(sample_641, fast=FALSE)
sample_641$high.tss <- ifelse(sample_641$TSS.enrichment > 2, 'High', 'Low')

VlnPlot(
  object = sample_641,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
sample_630 <- subset(
  x = sample_630,
  subset = nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
        nucleosome_signal < 2 &
    TSS.enrichment > 1
)
sample_630
#An object of class Seurat 
#149395 features across 5381 samples within 2 assays 
#Active assay: ATAC (117110 features, 0 variable features)
#1 other assay present: RNA

sample_636 <- subset(
  x = sample_636,
  subset = nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
sample_636
#An object of class Seurat 
#146215 features across 4241 samples within 2 assays 
#Active assay: ATAC (113930 features, 0 variable features)
#1 other assay present: RNA

sample_638 <- subset(
  x = sample_638,
  subset = nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
sample_638
#An object of class Seurat 
#152665 features across 4887 samples within 2 assays 
#Active assay: ATAC (120380 features, 0 variable features)
#1 other assay present: RNA

sample_641 <- subset(
  x = sample_641,
  subset = nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
sample_641
#An object of class Seurat 
#149036 features across 5328 samples within 2 assays 
#Active assay: ATAC (116751 features, 0 variable features)
#1 other assay present: RNA

#add metadata for each sample
sample_630$groupid <- 'adlib'
sample_636$groupid <- 'dehydrated'
sample_638$groupid <- 'dehydrated'
sample_641$groupid <- 'adlib'
sample_630$sex <- 'female'
sample_636$sex <- 'female'
sample_638$sex <- 'male'
sample_641$sex <- 'male'

#signac says we need to recall peaks with MACS2.  I have MACS2 loaded in my NGS environment
peaks.630 <- CallPeaks(sample_630, macs2.path = "/Users/kellyhyndman/NGS_env/bin/macs3")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks.630 <- keepStandardChromosomes(peaks.630, pruning.mode = "coarse")
peaks.630 <- subsetByOverlaps(x = peaks.630, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts.630 <- FeatureMatrix(
  fragments = Fragments(sample_630),
  features = peaks.630,
  cells = colnames(sample_630)
)

# create a new assay using the MACS3 peak set and add it to the Seurat object. the Assay Peaks is the new called peaks that are better than cellranger.
sample_630[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.630,
  fragments = frags.630,
  annotation = annotations,
)
sample_630

#An object of class Seurat 
#269495 features across 5381 samples within 3 assays 
#Active assay: ATAC (117110 features, 0 variable features)
#2 other assays present: RNA, peaks

peaks.636 <- CallPeaks(sample_636, macs2.path = "/Users/kellyhyndman/NGS_env/bin/macs3")
peaks.636 <- keepStandardChromosomes(peaks.636, pruning.mode = "coarse")
peaks.636 <- subsetByOverlaps(x = peaks.636, ranges = blacklist_mm10, invert = TRUE)
macs2_counts.636 <- FeatureMatrix(
  fragments = Fragments(sample_636),
  features = peaks.636,
  cells = colnames(sample_636)
)

sample_636[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.636,
  fragments = frags.636,
  annotation = annotations
)
sample_636
#An object of class Seurat 
#264740 features across 4241 samples within 3 assays 
#Active assay: ATAC (113930 features, 0 variable features)
#2 other assays present: RNA, peaks

peaks.638 <- CallPeaks(sample_638, macs2.path = "/Users/kellyhyndman/NGS_env/bin/macs3")
peaks.638 <- keepStandardChromosomes(peaks.638, pruning.mode = "coarse")
peaks.638 <- subsetByOverlaps(x = peaks.638, ranges = blacklist_mm10, invert = TRUE)
macs2_counts.638 <- FeatureMatrix(
  fragments = Fragments(sample_638),
  features = peaks.638,
  cells = colnames(sample_638)
)

sample_638[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.638,
  fragments = frags.638,
  annotation = annotations
)
sample_638

#An object of class Seurat 
#279315 features across 4887 samples within 3 assays 
#Active assay: ATAC (120380 features, 0 variable features)
#2 other assays present: RNA, peaks

peaks.641 <- CallPeaks(sample_641, macs2.path = "/Users/kellyhyndman/NGS_env/bin/macs3")
peaks.641 <- keepStandardChromosomes(peaks.641, pruning.mode = "coarse")
peaks.641 <- subsetByOverlaps(x = peaks.641, ranges = blacklist_mm10, invert = TRUE)
macs2_counts.641 <- FeatureMatrix(
  fragments = Fragments(sample_641),
  features = peaks.641,
  cells = colnames(sample_641)
)

sample_641[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.641,
  fragments = frags.641,
  annotation = annotations
)
sample_641

#An object of class Seurat 
#271318 features across 5328 samples within 3 assays 
#Active assay: ATAC (116751 features, 0 variable features)
#2 other assays present: RNA, peaks

DefaultAssay(sample_630) <- "RNA"
sample_630<-SCTransform(sample_630)
sample_638<-SCTransform(sample_638)
sample_636<-SCTransform(sample_636)
sample_641<-SCTransform(sample_641)





sample_630 <-NormalizeData(sample_630, normalization.method = "LogNormalize", scale.factor = 10000)
sample_630 <- FindVariableFeatures(sample_630, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample_630)
sample_630 <- ScaleData(sample_630, features = all.genes)
sample_630 <- RunPCA(sample_630, features = VariableFeatures(object = sample_630))
DimHeatmap(sample_630, dims = 1:15, cells = 500, balanced = TRUE)
sample_630 <- JackStraw(sample_630, num.replicate = 100)
sample_630 <- ScoreJackStraw(sample_630, dims = 1:20)
JackStrawPlot(sample_630, dims = 1:20)
sample_630 <- FindNeighbors(sample_630, dims = 1:20)
sample_630 <- FindClusters(sample_630, resolution = 1.4)
sample_630 <- RunUMAP(sample_630, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p630 <- DimPlot(sample_630, label = TRUE)
p630

DefaultAssay(sample_630) <- "peaks"
sample_630 <- RunTFIDF(sample_630)
sample_630 <- FindTopFeatures(sample_630, min.cutoff = 'q0')
sample_630 <- RunSVD(object = sample_630)
gene.activities <- GeneActivity(sample_630)
sample_630[['Gene']] <- CreateAssayObject(counts = gene.activities)
sample_630

DefaultAssay(sample_636) <- "RNA"
sample_636 <-NormalizeData(sample_636, normalization.method = "LogNormalize", scale.factor = 10000)
sample_636 <- FindVariableFeatures(sample_636, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample_636)
sample_636 <- ScaleData(sample_636, features = all.genes)
sample_636 <- RunPCA(sample_636, features = VariableFeatures(object = sample_636))
DimHeatmap(sample_636, dims = 1:15, cells = 500, balanced = TRUE)
sample_636 <- JackStraw(sample_636, num.replicate = 100)
sample_636 <- ScoreJackStraw(sample_636, dims = 1:20)
JackStrawPlot(sample_636, dims = 1:20)
sample_636 <- FindNeighbors(sample_636, dims = 1:20)
sample_636 <- FindClusters(sample_636, resolution = 1.4)
sample_636 <- RunUMAP(sample_636, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p636 <- DimPlot(sample_636, label = TRUE)
p636

DefaultAssay(sample_636) <- "peaks"
sample_636 <- RunTFIDF(sample_636)
sample_636 <- FindTopFeatures(sample_636, min.cutoff = 'q0')
sample_636 <- RunSVD(object = sample_636)
gene.activities <- GeneActivity(sample_636)
sample_636[['Gene']] <- CreateAssayObject(counts = gene.activities)
sample_636

DefaultAssay(sample_638) <- "RNA"
sample_638 <-NormalizeData(sample_638, normalization.method = "LogNormalize", scale.factor = 10000)
sample_638 <- FindVariableFeatures(sample_638, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample_638)
sample_638 <- ScaleData(sample_638, features = all.genes)
sample_638 <- RunPCA(sample_638, features = VariableFeatures(object = sample_638))
DimHeatmap(sample_638, dims = 1:15, cells = 500, balanced = TRUE)
sample_638 <- JackStraw(sample_638, num.replicate = 100)
sample_638 <- ScoreJackStraw(sample_638, dims = 1:20)
JackStrawPlot(sample_638, dims = 1:20)
sample_638 <- FindNeighbors(sample_638, dims = 1:20)
sample_638 <- FindClusters(sample_638, resolution = 1.4)
sample_638 <- RunUMAP(sample_638, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p638 <- DimPlot(sample_638, label = TRUE)
p638

DefaultAssay(sample_638) <- "peaks"
sample_638 <- RunTFIDF(sample_638)
sample_638 <- FindTopFeatures(sample_638, min.cutoff = 'q0')
sample_638 <- RunSVD(object = sample_638)
gene.activities <- GeneActivity(sample_638)
sample_638[['Gene']] <- CreateAssayObject(counts = gene.activities)
sample_638

DefaultAssay(sample_641) <- "RNA"
sample_641 <-NormalizeData(sample_641, normalization.method = "LogNormalize", scale.factor = 10000)
sample_641 <- FindVariableFeatures(sample_641, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample_641)
sample_641 <- ScaleData(sample_641, features = all.genes)
sample_641 <- RunPCA(sample_641, features = VariableFeatures(object = sample_641))
DimHeatmap(sample_641, dims = 1:15, cells = 500, balanced = TRUE)
sample_641 <- JackStraw(sample_641, num.replicate = 100)
sample_641 <- ScoreJackStraw(sample_641, dims = 1:20)
JackStrawPlot(sample_641, dims = 1:20)
sample_641 <- FindNeighbors(sample_641, dims = 1:20)
sample_641 <- FindClusters(sample_641, resolution = 1.4)
sample_641 <- RunUMAP(sample_641, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p641 <- DimPlot(sample_641, label = TRUE)
p641

DefaultAssay(sample_641) <- "peaks"
sample_641 <- RunTFIDF(sample_641)
sample_641 <- FindTopFeatures(sample_641, min.cutoff = 'q0')
sample_641 <- RunSVD(object = sample_641)
gene.activities <- GeneActivity(sample_641)
sample_641[['Gene']] <- CreateAssayObject(counts = gene.activities)
sample_641

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = sample_630,
  y = list(sample_636, sample_638, sample_641))

combined[["peaks"]]
#ChromatinAssay data with 152375 features for 19837 cells
#Variable features: 0 
#Genome: 
 # Annotation present: TRUE 
#Motifs present: FALSE 
#Fragment files: 4 

combined[["RNA"]]
#Assay data with 32285 features for 19837 cells
#First 10 features:
 # Xkr4, Gm1992, Gm19938, Gm37381, Rp1, Sox17, Gm37587, Gm37323, Mrpl15, Lypla1 

combined[["ATAC"]]
#ChromatinAssay data with 147370 features for 19837 cells
#Variable features: 0 
#Genome: 
#Annotation present: TRUE 
#Motifs present: FALSE 
#Fragment files: 4 

saveRDS(combined, file = "~/Desktop/Dehydration/combined.rds")
DimPlot(combined)
#Gene expression data processing
#We can normalize the gene expression data using SCTransform or lognorm, and reduce the dimensionality using PCA.

DefaultAssay(combined) <- "RNA"
combined <-NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE)
combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)
JackStrawPlot(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 1.4)
combined <- RunUMAP(combined, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p2 <- DimPlot(combined, label = TRUE, reduction = 'umap.rna')
p2
p3 <- DimPlot(combined, label = TRUE, reduction = 'lsi')
p3
combined.marker.RNA <- FindAllMarkers(combined, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.marker.RNA %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(combined.marker.RNA, file = "~/Desktop/Dehydration/RNAmarkers.csv")
saveRDS(combined, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/combined.rds")

#SCTtransform
DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 1.4)
combined <- RunUMAP(combined, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p1a <- DimPlot(combined, label = TRUE)
p1a
sctcombined.marker <- FindAllMarkers(combined, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sctcombined.marker, file = "~/Desktop/Dehydration/RNAmarkers.csv")
saveRDS(combined, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/combined.rds")

#Peaks normalization and scaling.
DefaultAssay(combined) <- "peaks"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(object = combined)
gene.activities <- GeneActivity(combined)
combined[['Gene']] <- CreateAssayObject(counts = gene.activities)

combined <- NormalizeData(
  object = combined,
  assay = 'Gene',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_Gene))
saveRDS(combined, file = "~/Desktop/Dehydration/combined.rds")


VariableFeatures(combined[["SCT"]]) <- rownames(combined[["SCT"]]@scale.data)
features <- SelectIntegrationFeatures(object.list = combined)


DefaultAssay(combined) <- "RNA"
DefaultAssay(combined) <- "peaks"
DefaultAssay(combined) <- "ATAC"
DefaultAssay(combined) <- "SCT"
DefaultAssay(combined) <- "Gene"
DimPlot(object = combined, label = TRUE, repel = TRUE)

DimPlot(combined, label = TRUE, group.by = "sex")


#DAC for clusters
DefaultAssay(dehydration) <- "peaks"
Idents(dehydration) <- "celltype"
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding dac for: ",cluster))
  dehydration <- seurat_aggregate
  dac <- FindMarkers(dehydration, 
                     ident.1 = cluster,  
                     logfc.threshold = 0.25,
                     test.use = 'LR', 
                     min.pct = 0.01,
                     only.pos = FALSE,
                     latent.vars = 'atac_peak_region_fragments')
  open <- rownames(dac)  
  cf <- ClosestFeature(dehydration, regions = open)
  return(cbind(dac, gene=cf$gene_name, distance=cf$distance))}
idents <- levels(dehydration@meta.data$celltype)
list.cluster.dac <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = dehydration)})
write.xlsx(list.cluster.dac, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/clusterdac.xlsx", sheetName = idents, rowNames = T)

write.csv(peakscombined.marker, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/peaksmarkers22.csv")
saveRDS(dehydration, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/dehydration.rds")

dac <- FindAllMarkers(dehydration, 
                   test.use = 'LR', 
                   min.pct = 0.1,
                   latent.vars = 'atac_peak_region_fragments')

write.csv(dac, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/peaksmarkers22.csv")



#ran files on local mac with Azimuth to get ideas of cell ideas by cross referencing the findallmarkers with azimuth.
#export a smaller file for Azimuth
DefaultAssay(dehydration) <- "RNA"
Azi <- DietSeurat(object = dehydration, assays = "RNA")
saveRDS(Azi, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/combineazi.rds")

#export this file and run Azimuth on webinterface.  Then manually cross reference with Cell cluster DEGs.
#now will rename idents
dehydration <- RenameIdents(dehydration, '0' = "fPTS3", '1' = "mPTS3", '2' = "fPTS1", '3' = "cTAL2", '4' = "Peritubular_EC", '5' = "fdTL", '6' = "mPTS2", '7' = "mPTS3", '8' = "mTAL", '9' = "mdTL", '10' = "fPTS2", '11' = "DCT1", '12' = "mPTS1",  '13' = "cTAL1", '14' = "DCT2", '15' = "mPTS1", '16' = "CDPC", '17' = "Fibroblast", '18' = "CNT2", '19' = "fPTS1", '20' = "VR1", '21' = "ICa", '22' = "dTL", '23'="cTAL3", '24'= "ICb", '25'='Monocytes', '26'="CNT1", '27' = "med_Fibro", '28' = "IMCD", "29"= "DCT2b", '30' = "CNT-PC", '31' = "Podocyte", '32' = "Tcell", '33'="Parietal", '34' = "VR2", '35' = "CNT3", '36'="VR3")
dehydration$celltype <- Idents(dehydration)
dehydration <- RenameIdents(dehydration, 'cTAL' = "cTAL2")
levels(dehydration) <- c("Podocyte", "Parietal", "mPTS1", "mPTS2", "mPTS3","fPTS1","fPTS2", "fPTS3", "mdTL", "fdTL", "dTL", "mTAL", "cTAL1","cTAL2", "cTAL3", "DCT1", "DCT2", "DCT2b", "CNT1", "CNT2", "CNT3", "CNT-PC", "CDPC", "IMCD", "ICa", "ICb", "Peritubular_EC", "VR1", "VR2", "VR3", "Fibroblast", "med_Fibro", "Monocytes", "Tcell")

p1 <-DimPlot(object = dehydration, label = TRUE, repel = TRUE)
p1

saveRDS(dehydration, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/dehydration.rds")

#get cell counts etc
table(Idents(dehydration))
table(dehydration$groupid)
table(dehydration$sex)
table(Idents(dehydration), dehydration$groupid)
table(Idents(dehydration), dehydration$sex)
prop.table(table(Idents(dehydration)))
prop.table(table(Idents(dehydration), dehydration$groupid), margin = 2)

#to change idents to separate con and KO. Original Clusters are in meta$seurat_clusters.  Then we can separate based on sex.
dehydration$celltype.groupid <- paste(Idents(dehydration), dehydration$groupid, sep = "_")
Idents(dehydration) <- "celltype.groupid"
dehydration$celltype.sex <- paste(Idents(dehydration), dehydration$sex, sep = "_")
Idents(dehydration) <- "celltype"
Idents(dehydration) <- "celltype.sex"
Idents(dehydration) <-"seurat_clusters"

DimPlot(object = dehydration, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = dehydration, label = TRUE, repel = TRUE) 

saveRDS(dehydration, file = "/data/user/hyndmank/multiomic/dehydration2021/multiome_Final/dehydration.rds")

#export barcode ids
Idents(dehydration) <- "celltype"
dehydrate_barcodes <-Idents(dehydration)
write.csv(dehydrate_barcodes, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/dehydratebarcodes.csv")

#DEG
DefaultAssay(dehydration) <- "RNA"
output <- "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/test/DES_"
Idents(dehydration) <- "celltype"
for (i in (levels(dehydration)))({
  try({
    ident1 <- paste0(i,"_adlib")
    ident2 <- paste0(i,"_dehydrated")
    Idents(dehydration) <- "celltype.groupid"
    deg <- FindMarkers(dehydration, ident.1 = ident1, ident.2=ident2, test.use = "DESeq2", verbose = TRUE) 
    write.csv(deg, file=paste0(output,i,".csv"))
  })
})


#DEG wilcox
DefaultAssay(dehydration) <- "RNA"
output <- "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DEG_Deseq/"
Idents(dehydration) <- "celltype"
for (i in (levels(dehydration))){
  try({
    ident1 <- paste0(i,"_adlib")
    ident2 <- paste0(i,"_dehydrated")
    Idents(dehydration) <- "celltype.groupid"
    degwilcox <- FindMarkers(dehydration, ident.1 = ident1, ident.2=ident2, verbose = TRUE) 
    write.csv(degwilcox, file=paste0(output,i,".csv"))
     })
}
#DAC had some issues with crashing so saving only the peaks to try running locally.
DefaultAssay(dehydration) <- "peaks"
peaks<-DietSeurat(
  dehydration,
  counts = TRUE,
  data = TRUE,
  scale.data = TRUE,
  features = NULL,
  assays = 'peaks',
  dimreducs = 'lsi',
  graphs = NULL
)
saveRDS(peaks, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/peaks.rds")

#DAC
DefaultAssay(dehydration) <- "peaks"
output <- "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_"
Idents(dehydration) <- "celltype"
for (i in (levels(dehydration))){
  try({
    ident1 <- paste0(i,"_adlib")
    ident2 <- paste0(i,"_dehydrated")
    Idents(dehydration) <- "celltype.groupid"
    DAC <- FindMarkers(dehydration, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(dehydration, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}

#DAC had some issues with crashing so saving only the peaks to try running locally.
DefaultAssay(dehydration) <- "peaks"
peaks<-DietSeurat(
  dehydration,
  counts = TRUE,
  data = TRUE,
  scale.data = TRUE,
  features = NULL,
  assays = 'peaks',
  dimreducs = 'lsi',
  graphs = NULL
)
saveRDS(peaks, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/peaks.rds")

#DAC females
DefaultAssay(peaks) <- "peaks"
setwd("~/Desktop")
output <- "~/Desktop/DAC_"
Idents(peaks) <- "celltype"
for (i in (levels(peaks))){
  try({
    ident1 <- paste0(i,"_adlib_female")
    ident2 <- paste0(i,"_dehydrated_female")
    Idents(peaks) <- "celltype.sex"
    DAC <- FindMarkers(peaks, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(peaks, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
setwd("~/Desktop/Fem_DAC")
list_of_files <- list.files(path = "~/Desktop/Fem_DAC/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <- subset(df, data$p_val_adj<0.05) 
new_df <-df %>% filter(p_val_adj <0.05)
write.xlsx(new_df, file = "~/Desktop/Fem_DAC/Fem_DAC.xlsx")

#DAC males
DefaultAssay(peaks) <- "peaks"
setwd("/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male")
output <- "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/M_"
Idents(peaks) <- "celltype"
for (i in (levels(peaks))){
  try({
    ident1 <- paste0(i,"_adlib_male")
    ident2 <- paste0(i,"_dehydrated_male")
    Idents(peaks) <- "celltype.sex"
    DAC <- FindMarkers(peaks, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(peaks, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <- subset(df, data$p_val_adj<0.05) 
new_df <-df %>% filter(p_val_adj <0.05)
write.xlsx(new_df, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_DAC.xlsx")





#########################################################################################################


DefaultAssay(dehydration) <- "peaks"
Idents(dehydration) <- "celltype.groupid"

DiffGenesCDPC <- FindMarkers(dehydration, ident.1 = "CDPC_adlib", ident.2 = "CDPC_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCDPC)  
cf <- ClosestFeature(dehydration, regions = open)
CDPC <-(cbind(DiffGenesCDPC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CDPC, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_CDPC.xlsx")

DiffGenesCNTPC <- FindMarkers(dehydration, ident.1 = "CNT-PC_adlib", ident.2 = "CNT-PC_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNTPC)  
cf <- ClosestFeature(dehydration, regions = open)
CNTPC <-(cbind(DiffGenesCNTPC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNTPC, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_CNTPC.xlsx")

DiffGenesCNT1 <- FindMarkers(dehydration, ident.1 = "CNT1_adlib", ident.2 = "CNT1_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT1)  
cf <- ClosestFeature(dehydration, regions = open)
CNT1 <-(cbind(DiffGenesCNT1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_CNT1.xlsx")

DiffGenesCNT2 <- FindMarkers(dehydration, ident.1 = "CNT2_adlib", ident.2 = "CNT2_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT2)  
cf <- ClosestFeature(dehydration, regions = open)
CNT2 <-(cbind(DiffGenesCNT2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_CNT2.xlsx")

DiffGenesCNT3 <- FindMarkers(dehydration, ident.1 = "CNT3_adlib", ident.2 = "CNT3_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT3)  
cf <- ClosestFeature(dehydration, regions = open)
CNT3 <-(cbind(DiffGenesCNT3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_CNT3.xlsx")

DiffGenesICb <- FindMarkers(dehydration, ident.1 = "ICb_adlib", ident.2 = "ICb_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesICb)  
cf <- ClosestFeature(dehydration, regions = open)
ICb <-(cbind(DiffGenesICb, gene=cf$gene_name, distance=cf$distance))
write.xlsx(ICb, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_ICb.xlsx")

DiffGenesICa <- FindMarkers(dehydration, ident.1 = "ICa_adlib", ident.2 = "ICa_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesICa)  
cf <- ClosestFeature(dehydration, regions = open)
ICa <-(cbind(DiffGenesICa, gene=cf$gene_name, distance=cf$distance))
write.xlsx(ICa, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_ICa.xlsx")

DiffGenesIMCD <- FindMarkers(dehydration, ident.1 = "IMCD_adlib", ident.2 = "IMCD_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesIMCD)  
cf <- ClosestFeature(dehydration, regions = open)
IMCD <-(cbind(DiffGenesIMCD, gene=cf$gene_name, distance=cf$distance))
write.xlsx(IMCD, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_IMCD.xlsx")

DiffGenesFibroblast <- FindMarkers(dehydration, ident.1 = "Fibroblast_adlib", ident.2 = "Fibroblast_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesFibroblast)  
cf <- ClosestFeature(dehydration, regions = open)
Fibroblast <-(cbind(DiffGenesFibroblast, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Fibroblast, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_Fibroblast.xlsx")

DiffGenesPeritubular_EC <- FindMarkers(dehydration, ident.1 = "Peritubular_EC_adlib", ident.2 = "Peritubular_EC_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesPeritubular_EC)  
cf <- ClosestFeature(dehydration, regions = open)
Peritubular_EC <-(cbind(DiffGenesPeritubular_EC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Peritubular_EC, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_Peritubular_EC.xlsx")

DiffGenesVR1 <- FindMarkers(dehydration, ident.1 = "VR1_adlib", ident.2 = "VR1_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesVR1)  
cf <- ClosestFeature(dehydration, regions = open)
VR1 <-(cbind(DiffGenesVR1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(VR1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_VR1.xlsx")

DiffGenesmed_Fibro <- FindMarkers(dehydration, ident.1 = "med_Fibro_adlib", ident.2 = "med_Fibro_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmed_Fibro)  
cf <- ClosestFeature(dehydration, regions = open)
med_Fibro <-(cbind(DiffGenesmed_Fibro, gene=cf$gene_name, distance=cf$distance))
write.xlsx(med_Fibro, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_med_Fibro.xlsx")

DiffGenesMonocytes <- FindMarkers(dehydration, ident.1 = "Monocytes_adlib", ident.2 = "Monocytes_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesMonocytes)  
cf <- ClosestFeature(dehydration, regions = open)
Monocytes <-(cbind(DiffGenesMonocytes, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Monocytes, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_Monocytes.xlsx")

DiffGenesTcell <- FindMarkers(dehydration, ident.1 = "Tcell_adlib", ident.2 = "Tcell_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesTcell)  
cf <- ClosestFeature(dehydration, regions = open)
Tcell <-(cbind(DiffGenesTcell, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Tcell, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_Tcell.xlsx")


DiffGenesPOD <- FindMarkers(dehydration, ident.1 = "Podocyte_adlib", ident.2 = "Podocyte_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesPOD)  
cf <- ClosestFeature(dehydration, regions = open)
Podocyte <-(cbind(DiffGenesPOD, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Podocyte, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_Podocyte.xlsx")


DiffGenesfPTS1 <- FindMarkers(dehydration, ident.1 = "fPTS1_adlib", ident.2 = "fPTS1_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfPTS1)  
cf <- ClosestFeature(dehydration, regions = open)
fPTS1 <-(cbind(DiffGenesfPTS1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_fPTS1.xlsx")

DiffGenesmPTS1 <- FindMarkers(dehydration, ident.1 = "mPTS1_adlib", ident.2 = "mPTS1_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmPTS1)  
cf <- ClosestFeature(dehydration, regions = open)
mPTS1 <-(cbind(DiffGenesmPTS1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_mPTS1.xlsx")

DiffGenesmPTS2 <- FindMarkers(dehydration, ident.1 = "mPTS2_adlib", ident.2 = "mPTS2_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmPTS2)  
cf <- ClosestFeature(dehydration, regions = open)
mPTS2 <-(cbind(DiffGenesmPTS2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_mPTS2.xlsx")

DiffGenesfPTS2 <- FindMarkers(dehydration, ident.1 = "fPTS2_adlib", ident.2 = "fPTS2_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfPTS2)  
cf <- ClosestFeature(dehydration, regions = open)
fPTS2 <-(cbind(DiffGenesfPTS2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_fPTS2.xlsx")

DiffGenesfPTS3 <- FindMarkers(dehydration, ident.1 = "fPTS3_adlib", ident.2 = "fPTS3_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfPTS3)  
cf <- ClosestFeature(dehydration, regions = open)
fPTS3 <-(cbind(DiffGenesfPTS3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_fPTS3.xlsx")

DiffGenesmPTS3 <- FindMarkers(dehydration, ident.1 = "mPTS3_adlib", ident.2 = "mPTS3_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmPTS3)  
cf <- ClosestFeature(dehydration, regions = open)
mPTS3 <-(cbind(DiffGenesmPTS3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_mPTS3.xlsx")

DiffGenesfdTL <- FindMarkers(dehydration, ident.1 = "fdTL_adlib", ident.2 = "fdTL_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfdTL)  
cf <- ClosestFeature(dehydration, regions = open)
fdTL <-(cbind(DiffGenesfdTL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fdTL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_fdTL.xlsx")

DiffGenesmdTL <- FindMarkers(dehydration, ident.1 = "mdTL_adlib", ident.2 = "mdTL_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmdTL)  
cf <- ClosestFeature(dehydration, regions = open)
mdTL <-(cbind(DiffGenesmdTL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mdTL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_mdTL.xlsx")

DiffGenesdTL <- FindMarkers(dehydration, ident.1 = "dTL_adlib", ident.2 = "dTL_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesdTL)  
cf <- ClosestFeature(dehydration, regions = open)
dTL <-(cbind(DiffGenesdTL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(dTL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_dTL.xlsx")

DiffGenesmTAL <- FindMarkers(dehydration, ident.1 = "mTAL_adlib", ident.2 = "mTAL_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmTAL)  
cf <- ClosestFeature(dehydration, regions = open)
mTAL <-(cbind(DiffGenesmTAL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mTAL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_mTAL.xlsx")

DiffGenescTAL2 <- FindMarkers(dehydration, ident.1 = "cTAL2_adlib", ident.2 = "cTAL2_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenescTAL2)  
cf <- ClosestFeature(dehydration, regions = open)
cTAL <-(cbind(DiffGenescTAL2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_cTAL2.xlsx")

DiffGenescTAL1 <- FindMarkers(dehydration, ident.1 = "cTAL1_adlib", ident.2 = "cTAL1_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenescTAL1)  
cf <- ClosestFeature(dehydration, regions = open)
cTAL1 <-(cbind(DiffGenescTAL1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_cTAL1.xlsx")

DiffGenescTAL3 <- FindMarkers(dehydration, ident.1 = "cTAL3_adlib", ident.2 = "cTAL3_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenescTAL3)  
cf <- ClosestFeature(dehydration, regions = open)
cTAL3 <-(cbind(DiffGenescTAL3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_cTAL3.xlsx")

DiffGenesDCT1 <- FindMarkers(dehydration, ident.1 = "DCT1_adlib", ident.2 = "DCT1_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesDCT1)  
cf <- ClosestFeature(dehydration, regions = open)
DCT1 <-(cbind(DiffGenesDCT1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_DCT1.xlsx")

DiffGenesDCT2 <- FindMarkers(dehydration, ident.1 = "DCT2_adlib", ident.2 = "DCT2_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesDCT2)  
cf <- ClosestFeature(dehydration, regions = open)
DCT2 <-(cbind(DiffGenesDCT2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_DCT2.xlsx")

DiffGenesDCT2b <- FindMarkers(dehydration, ident.1 = "DCT2b_adlib", ident.2 = "DCT2b_dehydrated", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesDCT2b)  
cf <- ClosestFeature(dehydration, regions = open)
DCT2b <-(cbind(DiffGenesDCT2b, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT2b, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_DCT2b.xlsx")

#problem ones
DiffGenesVR2 <- FindMarkers(dehydration, ident.1 = "VR2_adlib", ident.2 = "VR2_dehydrated", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesVR2)  
cf <- ClosestFeature(dehydration, regions = open)
VR2 <-(cbind(DiffGenesVR2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(VR2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_VR2.xlsx")

DiffGenesVR3 <- FindMarkers(dehydration, ident.1 = "VR3_adlib", ident.2 = "VR3_dehydrated", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesVR3)  
cf <- ClosestFeature(dehydration, regions = open)
VR3 <-(cbind(DiffGenesVR3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(VR3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_VR3.xlsx")

DiffGenesParietal <- FindMarkers(dehydration, ident.1 = "Parietal_adlib", ident.2 = "Parietal_dehydrated", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesParietal)  
cf <- ClosestFeature(dehydration, regions = open)
Parietal <-(cbind(DiffGenesParietal, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Parietal, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_Parietal.xlsx")

DiffGenesCNT3 <- FindMarkers(dehydration, ident.1 = "CNT3_adlib", ident.2 = "CNT3_dehydrated", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT3)  
cf <- ClosestFeature(dehydration, regions = open)
CNT3 <-(cbind(DiffGenesCNT3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/DAC_CNT3.xlsx")


#DEGfor sex specific

DefaultAssay(dehydration) <- "RNA"
output <- "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/sex/"
Idents(dehydration) <- "celltype"
for (i in (levels(dehydration)))({
  try({
    ident1 <- paste0(i,"_adlib_male")
    ident2 <- paste0(i,"_dehydrated_male")
    Idents(dehydration) <- "celltype.sex"
    deg <- FindMarkers(dehydration, ident.1 = ident1, ident.2=ident2, test.use = "DESeq2", verbose = TRUE) 
    write.csv(deg, file=paste0(output,i,".csv"))
  })
})

DefaultAssay(dehydration) <- "RNA"
output <- "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/sex/female/"
Idents(dehydration) <- "celltype"
for (i in (levels(dehydration)))({
  try({
    ident1 <- paste0(i,"_adlib_female")
    ident2 <- paste0(i,"_dehydrated_female")
    Idents(dehydration) <- "celltype.sex"
    deg <- FindMarkers(dehydration, ident.1 = ident1, ident.2=ident2, test.use = "DESeq2", verbose = TRUE) 
    write.csv(deg, file=paste0(output,i,".csv"))
  })
})

#DAC per sex.  full script crashed so ran individually
DefaultAssay(dehydration) <- "peaks"
Idents(dehydration) <- "celltype.sex"

DiffGenesCDPC <- FindMarkers(dehydration, ident.1 = "CDPC_adlib_male", ident.2 = "CDPC_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCDPC)  
cf <- ClosestFeature(dehydration, regions = open)
CDPC <-(cbind(DiffGenesCDPC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CDPC, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_CDPC.xlsx")


DiffGenesCNT1 <- FindMarkers(dehydration, ident.1 = "CNT1_adlib_male", ident.2 = "CNT1_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT1)  
cf <- ClosestFeature(dehydration, regions = open)
CNT1 <-(cbind(DiffGenesCNT1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_CNT1.xlsx")

DiffGenesCNT2 <- FindMarkers(dehydration, ident.1 = "CNT2_adlib_male", ident.2 = "CNT2_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT2)  
cf <- ClosestFeature(dehydration, regions = open)
CNT2 <-(cbind(DiffGenesCNT2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_CNT2.xlsx")

DiffGenesCNT3 <- FindMarkers(dehydration, ident.1 = "CNT3_adlib_male", ident.2 = "CNT3_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT3)  
cf <- ClosestFeature(dehydration, regions = open)
CNT3 <-(cbind(DiffGenesCNT3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_CNT3.xlsx")

DiffGenesICb <- FindMarkers(dehydration, ident.1 = "ICb_adlib_male", ident.2 = "ICb_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesICb)  
cf <- ClosestFeature(dehydration, regions = open)
ICb <-(cbind(DiffGenesICb, gene=cf$gene_name, distance=cf$distance))
write.xlsx(ICb, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_ICb.xlsx")

DiffGenesICa <- FindMarkers(dehydration, ident.1 = "ICa_adlib_male", ident.2 = "ICa_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesICa)  
cf <- ClosestFeature(dehydration, regions = open)
ICa <-(cbind(DiffGenesICa, gene=cf$gene_name, distance=cf$distance))
write.xlsx(ICa, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_ICa.xlsx")

DiffGenesIMCD <- FindMarkers(dehydration, ident.1 = "IMCD_adlib_male", ident.2 = "IMCD_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesIMCD)  
cf <- ClosestFeature(dehydration, regions = open)
IMCD <-(cbind(DiffGenesIMCD, gene=cf$gene_name, distance=cf$distance))
write.xlsx(IMCD, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_IMCD.xlsx")

DiffGenesFibroblast <- FindMarkers(dehydration, ident.1 = "Fibroblast_adlib_male", ident.2 = "Fibroblast_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesFibroblast)  
cf <- ClosestFeature(dehydration, regions = open)
Fibroblast <-(cbind(DiffGenesFibroblast, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Fibroblast, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_Fibroblast.xlsx")

DiffGenesPeritubular_EC <- FindMarkers(dehydration, ident.1 = "Peritubular_EC_adlib_male", ident.2 = "Peritubular_EC_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesPeritubular_EC)  
cf <- ClosestFeature(dehydration, regions = open)
Peritubular_EC <-(cbind(DiffGenesPeritubular_EC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Peritubular_EC, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_Peritubular_EC.xlsx")

DiffGenesVR1 <- FindMarkers(dehydration, ident.1 = "VR1_adlib_male", ident.2 = "VR1_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesVR1)  
cf <- ClosestFeature(dehydration, regions = open)
VR1 <-(cbind(DiffGenesVR1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(VR1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_VR1.xlsx")

DiffGenesmed_Fibro <- FindMarkers(dehydration, ident.1 = "med_Fibro_adlib_male", ident.2 = "med_Fibro_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmed_Fibro)  
cf <- ClosestFeature(dehydration, regions = open)
med_Fibro <-(cbind(DiffGenesmed_Fibro, gene=cf$gene_name, distance=cf$distance))
write.xlsx(med_Fibro, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_med_Fibro.xlsx")

DiffGenesMonocytes <- FindMarkers(dehydration, ident.1 = "Monocytes_adlib_male", ident.2 = "Monocytes_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesMonocytes)  
cf <- ClosestFeature(dehydration, regions = open)
Monocytes <-(cbind(DiffGenesMonocytes, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Monocytes, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_Monocytes.xlsx")

DiffGenesTcell <- FindMarkers(dehydration, ident.1 = "Tcell_adlib_male", ident.2 = "Tcell_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesTcell)  
cf <- ClosestFeature(dehydration, regions = open)
Tcell <-(cbind(DiffGenesTcell, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Tcell, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_Tcell.xlsx")

DiffGenesPOD <- FindMarkers(dehydration, ident.1 = "Podocyte_adlib_male", ident.2 = "Podocyte_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesPOD)  
cf <- ClosestFeature(dehydration, regions = open)
Podocyte <-(cbind(DiffGenesPOD, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Podocyte, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_Podocyte.xlsx")

DiffGenesfPTS1 <- FindMarkers(dehydration, ident.1 = "fPTS1_adlib_male", ident.2 = "fPTS1_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfPTS1)  
cf <- ClosestFeature(dehydration, regions = open)
fPTS1 <-(cbind(DiffGenesfPTS1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_fPTS1.xlsx")

DiffGenesmPTS1 <- FindMarkers(dehydration, ident.1 = "mPTS1_adlib_male", ident.2 = "mPTS1_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmPTS1)  
cf <- ClosestFeature(dehydration, regions = open)
mPTS1 <-(cbind(DiffGenesmPTS1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_mPTS1.xlsx")

DiffGenesmPTS2 <- FindMarkers(dehydration, ident.1 = "mPTS2_adlib_male", ident.2 = "mPTS2_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmPTS2)  
cf <- ClosestFeature(dehydration, regions = open)
mPTS2 <-(cbind(DiffGenesmPTS2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_mPTS2.xlsx")

DiffGenesfPTS2 <- FindMarkers(dehydration, ident.1 = "fPTS2_adlib_male", ident.2 = "fPTS2_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfPTS2)  
cf <- ClosestFeature(dehydration, regions = open)
fPTS2 <-(cbind(DiffGenesfPTS2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_fPTS2.xlsx")

DiffGenesfPTS3 <- FindMarkers(dehydration, ident.1 = "fPTS3_adlib_male", ident.2 = "fPTS3_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfPTS3)  
cf <- ClosestFeature(dehydration, regions = open)
fPTS3 <-(cbind(DiffGenesfPTS3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_fPTS3.xlsx")

DiffGenesmPTS3 <- FindMarkers(dehydration, ident.1 = "mPTS3_adlib_male", ident.2 = "mPTS3_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmPTS3)  
cf <- ClosestFeature(dehydration, regions = open)
mPTS3 <-(cbind(DiffGenesmPTS3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_mPTS3.xlsx")

DiffGenesfdTL <- FindMarkers(dehydration, ident.1 = "fdTL_adlib_male", ident.2 = "fdTL_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesfdTL)  
cf <- ClosestFeature(dehydration, regions = open)
fdTL <-(cbind(DiffGenesfdTL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fdTL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_fdTL.xlsx")

DiffGenesmdTL <- FindMarkers(dehydration, ident.1 = "mdTL_adlib_male", ident.2 = "mdTL_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmdTL)  
cf <- ClosestFeature(dehydration, regions = open)
mdTL <-(cbind(DiffGenesmdTL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mdTL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_mdTL.xlsx")

DiffGenesdTL <- FindMarkers(dehydration, ident.1 = "dTL_adlib_male", ident.2 = "dTL_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesdTL)  
cf <- ClosestFeature(dehydration, regions = open)
dTL <-(cbind(DiffGenesdTL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(dTL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_dTL.xlsx")

DiffGenesmTAL <- FindMarkers(dehydration, ident.1 = "mTAL_adlib_male", ident.2 = "mTAL_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesmTAL)  
cf <- ClosestFeature(dehydration, regions = open)
mTAL <-(cbind(DiffGenesmTAL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mTAL, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_mTAL.xlsx")

DiffGenescTAL2 <- FindMarkers(dehydration, ident.1 = "cTAL2_adlib_male", ident.2 = "cTAL2_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenescTAL2)  
cf <- ClosestFeature(dehydration, regions = open)
cTAL <-(cbind(DiffGenescTAL2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_cTAL2.xlsx")

DiffGenescTAL1 <- FindMarkers(dehydration, ident.1 = "cTAL1_adlib_male", ident.2 = "cTAL1_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenescTAL1)  
cf <- ClosestFeature(dehydration, regions = open)
cTAL1 <-(cbind(DiffGenescTAL1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_cTAL1.xlsx")

DiffGenescTAL3 <- FindMarkers(dehydration, ident.1 = "cTAL3_adlib_male", ident.2 = "cTAL3_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenescTAL3)  
cf <- ClosestFeature(dehydration, regions = open)
cTAL3 <-(cbind(DiffGenescTAL3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_cTAL3.xlsx")

DiffGenesDCT1 <- FindMarkers(dehydration, ident.1 = "DCT1_adlib_male", ident.2 = "DCT1_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesDCT1)  
cf <- ClosestFeature(dehydration, regions = open)
DCT1 <-(cbind(DiffGenesDCT1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT1, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_DCT1.xlsx")

DiffGenesDCT2 <- FindMarkers(dehydration, ident.1 = "DCT2_adlib_male", ident.2 = "DCT2_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesDCT2)  
cf <- ClosestFeature(dehydration, regions = open)
DCT2 <-(cbind(DiffGenesDCT2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_DCT2.xlsx")

DiffGenesDCT2b <- FindMarkers(dehydration, ident.1 = "DCT2b_adlib_male", ident.2 = "DCT2b_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesDCT2b)  
cf <- ClosestFeature(dehydration, regions = open)
DCT2b <-(cbind(DiffGenesDCT2b, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT2b, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_DCT2b.xlsx")

DiffGenesVR2 <- FindMarkers(dehydration, ident.1 = "VR2_adlib_male", ident.2 = "VR2_dehydrated_male", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesVR2)  
cf <- ClosestFeature(dehydration, regions = open)
VR2 <-(cbind(DiffGenesVR2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(VR2, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_VR2.xlsx")

DiffGenesVR3 <- FindMarkers(dehydration, ident.1 = "VR3_adlib_male", ident.2 = "VR3_dehydrated_male", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesVR3)  
cf <- ClosestFeature(dehydration, regions = open)
VR3 <-(cbind(DiffGenesVR3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(VR3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_VR3.xlsx")

DiffGenesParietal <- FindMarkers(dehydration, ident.1 = "Parietal_adlib_male", ident.2 = "Parietal_dehydrated_male", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesParietal)  
cf <- ClosestFeature(dehydration, regions = open)
Parietal <-(cbind(DiffGenesParietal, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Parietal, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_Parietal.xlsx")

DiffGenesCNT3 <- FindMarkers(dehydration, ident.1 = "CNT3_adlib_male", ident.2 = "CNT3_dehydrated_male", min.pct = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNT3)  
cf <- ClosestFeature(dehydration, regions = open)
CNT3 <-(cbind(DiffGenesCNT3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNT3, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_CNT3.xlsx")

DiffGenesCNTPC <- FindMarkers(dehydration, ident.1 = "CNT-PC_adlib_male", ident.2 = "CNT-PC_dehydrated_male", min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
open <- rownames(DiffGenesCNTPC)  
cf <- ClosestFeature(dehydration, regions = open)
CNTPC <-(cbind(DiffGenesCNTPC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CNTPC, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/DAC/male/male_CNTPC.xlsx")
####################################################
#DAC
DefaultAssay(dehydration) <- "peaks"
setwd("/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/test")
output <- ("/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/test/")
Idents(dehydration) <- "celltype"
for (i in (levels(dehydration))){
  try({
    ident1 <- paste0(i,"_adlib")
    ident2 <- paste0(i,"_dehydrated")
    Idents(dehydration) <- "celltype.groupid"
    DAC <- FindMarkers(dehydration, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0.05, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(dehydration, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
list_of_files <- list.files(path = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/test/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <-df %>% filter(p_val_adj <0.05)
write.xlsx(new_df, file = "/data/user/hyndmank/multiomic/dehydration2021/22_FINAL/test/All_DAC.xlsx")




DefaultAssay(peaks) <- "peaks"
setwd("~/Desktop")
output <- "~/Desktop/DAC_"
Idents(peaks) <- "celltype"
for (i in (levels(peaks))){
  try({
    ident1 <- paste0(i,"_adlib_female")
    ident2 <- paste0(i,"_dehydrated_female")
    Idents(peaks) <- "celltype.sex"
    DAC <- FindMarkers(peaks, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(peaks, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
setwd("~/Desktop/Fem_DAC")
list_of_files <- list.files(path = "~/Desktop/Fem_DAC/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <-df %>% filter(p_val_adj <0.05)
write.xlsx(new_df, file = "~/Desktop/Fem_DAC/Fem_DAC.xlsx")


DefaultAssay(peaks) <- "peaks"
output <- "~/Desktop/male_DAC/"
Idents(peaks) <- "celltype"
for (i in (levels(peaks))){
  try({
    ident1 <- paste0(i,"_adlib_male")
    ident2 <- paste0(i,"_dehydrated_male")
    Idents(peaks) <- "celltype.sex"
    DAC <- FindMarkers(peaks, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(peaks, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
setwd("~/Desktop/male_DAC")
list_of_files <- list.files(path = "~/Desktop/male_DAC/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <-df %>% filter(p_val_adj <0.05)
write.xlsx(new_df, file = "~/Desktop/male_DAC/male_DAC.xlsx")

write.xlsx(df, file = "~/Desktop/male_DAC/male_DAC_all.xlsx")

#####all figure code is in Figures.R in the Images folder in 22_FINAL
df <-data.frame(peaksmarkers22, stringsAsFactors = TRUE)
df2 <- df[,-1]
rownames(df2) <- df[,1]
open <-rownames(df2)

cf <- ClosestFeature(dehydration, regions = open)


#Tried the integration code, and it has good overlap, but poor # of clusters.  So not using it.
# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(combined, split.by = "groupid")
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "groupid")
p1 

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca",
                                         k.anchor = 20)
immune.combined <- IntegrateData(anchorset = immune.anchors)

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 1.4)
DimPlot(immune.combined, reduction = "umap")
