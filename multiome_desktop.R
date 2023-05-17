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

combined[["peaks"]]


#SCTtransform
DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 1.4)
combined <- RunUMAP(combined, dims = 1:20, reduction.name = 'umap.sct', reduction.key = 'rnaUMAP_')
p1a <- DimPlot(combined, label = TRUE)
p1a
sctcombined.marker <- FindAllMarkers(combined, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sctcombined.marker, file = "~/Desktop/Dehydration/RNAmarkers.csv")
saveRDS(combined, file = "~/Desktop/Dehydration/combined.rds")


#Gene expression data processing
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

p1b <- DimPlot(combined, label = TRUE, reduction = 'umap.rna')
p1b

DimPlot(combined, label = FALSE, reduction = 'umap.sct', group.by = 'sex')
DimPlot(combined, label = FALSE, reduction = 'umap.sct', group.by = 'groupid') 


combined.marker.RNA <- FindAllMarkers(combined, assay = "RNA", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
combined.marker.RNA %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(combined.marker.RNA, file = "~/Desktop/Dehydration/RNAmarkers2.csv")
saveRDS(combined, file = "~/Desktop/Dehydration/combined.rds")

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

DefaultAssay(combined) <- "RNA"
DefaultAssay(combined) <- "peaks"
DefaultAssay(combined) <- "ATAC"
DefaultAssay(combined) <- "SCT"
DefaultAssay(combined) <- "Gene"

#export this file and run Azimuth on webinterface.  Then manually cross reference with Cell cluster DEGs.
#now will rename idents
combined <- RenameIdents(combined, '0' = "fPTS3", '1' = "fPTS1", '2' = "cTAL2", '3' = "Peritubular_EC", '4' = "mPTS1", '5' = "mPTS2", '6' = "mdTL", '7' = "fPTS2", '8' = "mPTS3", '9' = "mTAL", '10' = "DCT1", '11' = "CNT1", '12' = "cTAL1",  '13' = "fdTL", '14' = "DCT2", '15' = "mPTS3", '16' = "mPTS3", '17' = "CDPC", '18' = "Fibroblast", '19' = "mPTS1", '20' = "fPTS1", '21' = "VR1", '22' = "ICa", '23'="dTL", '24'= "cTAL3", '25'='ICb', '26'="Monocytes", '27' = "fdTL", '28' = "med_Fibro", "29"= "IMCD", '30' = "DCT2b", '31' = "CNT-PC", '32' = "Podocyte", '33'="Tcell", '34' = "Parietal", '35' = "VR2", '36'="VR3", '37' = 'CNT2')
combined$celltype <- Idents(combined)

levels(combined) <- c("Podocyte", "Parietal", "mPTS1", "mPTS2", "mPTS3","fPTS1","fPTS2", "fPTS3", "mdTL", "fdTL", "dTL", "mTAL", "cTAL1","cTAL2", "cTAL3", "DCT1", "DCT2", "DCT2b", "CNT1", "CNT2", "CNT-PC", "CDPC", "IMCD", "ICa", "ICb", "Peritubular_EC", "VR1", "VR2", "VR3", "Fibroblast", "med_Fibro", "Monocytes", "Tcell")
combined$celltype <- Idents(combined)
DimPlot(object = combined, label = TRUE, repel = TRUE)
p1

saveRDS(combined, file = "~/Desktop/Dehydration/combined.rds")

#get cell counts etc
table(Idents(combined))
table(combined$groupid)
table(combined$sex)
table(Idents(combined), combined$groupid)
table(Idents(combined), combined$sex)
prop.table(table(Idents(combined)))
prop.table(table(Idents(combined), combined$groupid), margin = 2)

#to change idents to separate con and KO. Original Clusters are in meta$seurat_clusters.  Then we can separate based on sex.
combined$celltype.groupid <- paste(Idents(combined), combined$groupid, sep = "_")
Idents(combined) <- "celltype.groupid"
combined$celltype.sex <- paste(Idents(combined), combined$sex, sep = "_")
Idents(combined) <- "celltype"
Idents(combined) <- "celltype.sex"
combined$celltype.groupid.sex <- paste(Idents(combined), combined$sex, sep = "_")
Idents(combined) <- "celltype.groupid.sex"
Idents(combined) <-"seurat_clusters"

DimPlot(object = combined, label = TRUE, repel = TRUE, reduction = 'umap.sct') + NoLegend()
DimPlot(object = combined, label = TRUE, repel = TRUE) 



#DAC for clusters
DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype"
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding dac for: ",cluster))
  combined <- seurat_aggregate
  dac <- FindMarkers(combined, 
                     ident.1 = cluster,  
                     logfc.threshold = 0.25,
                     test.use = 'LR', 
                     min.pct = 0.1,
                     only.pos = FALSE,
                     latent.vars = 'atac_peak_region_fragments')
  open <- rownames(dac)  
  cf <- ClosestFeature(combined, regions = open)
  return(cbind(dac, gene=cf$gene_name, distance=cf$distance))}
idents <- levels(combined@meta.data$celltype)
list.cluster.dac <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = combined)})
write.xlsx(list.cluster.dac, file = "~/Desktop/Dehydration/clusterdac.xlsx", sheetName = idents, rowNames = T)

#to change excel tabs to a sheet with list
library(readxl)
library(tidyverse)

# specifying the path for file
path <- "~/Desktop/Dehydration/"

# set the working directory 
setwd(path)

# accessing all the sheets 
sheet = excel_sheets("clusterdac.xlsx")

# applying sheet names to dataframe names
data_frame = lapply(setNames(sheet, sheet), 
                    function(x) read_excel("clusterdac.xlsx", sheet=x))
# attaching all dataframes together
data_frame = bind_rows(data_frame, .id="Sheet")

# printing data of all sheets
print (data_frame)
write.xlsx(data_frame, file = "~/Desktop/Dehydration/clusterdac2.xlsx", rowNames = T)


#get cell counts etc
table(Idents(combined))
table(combined$groupid)
table(combined$sex)
table(Idents(combined), combined$groupid)
table(Idents(combined), combined$sex)
prop.table(table(Idents(combined)))
prop.table(table(Idents(combined), combined$groupid), margin = 2)

#to change idents to separate con and KO. Original Clusters are in meta$seurat_clusters.  Then we can separate based on sex.
combined$celltype.groupid <- paste(Idents(combined), combined$groupid, sep = "_")
Idents(combined) <- "celltype.groupid"
combined$celltype.sex <- paste(Idents(combined), combined$sex, sep = "_")
Idents(combined) <- "celltype"
Idents(combined) <- "celltype.sex"
Idents(combined) <-"seurat_clusters"

DimPlot(object = combined, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = combined, label = TRUE, repel = TRUE) 

saveRDS(combined, file = "/data/user/hyndmank/multiomic/dehydration2021/multiome_Final/dehydration.rds")

#export barcode ids
Idents(combined) <- "celltype"
dehydrate_barcodes <-Idents(combined)
write.csv(dehydrate_barcodes, file = "~/Desktop/combined/DEG/dehydratebarcodes.csv")

#DEG
DefaultAssay(combined) <- "RNA"
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined)))({
  try({
    ident1 <- paste0(i,"_dehydrated")
    ident2 <- paste0(i,"_adlib")
    Idents(combined) <- "celltype.groupid"
    deg <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
    deg <-cbind(deg, gene=rownames(deg))
    write.xlsx(deg,file=paste0(output,i,".xlsx"))
  })
})

list_of_files <- list.files(path = "~/Desktop/combined/DEG",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/DEG/dehy_DEG.xlsx")

#FEMALES
DefaultAssay(combined) <- "RNA"
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined)))({
  try({
    ident1 <- paste0(i,"_dehydrated_female")
    ident2 <- paste0(i,"_adlib_female")
    Idents(combined) <- "celltype.groupid.sex"
    deg <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
    deg <-cbind(deg, gene=rownames(deg))
    write.xlsx(deg,file=paste0(output,i,".xlsx"))
  })
})

list_of_files <- list.files(path = "~/Desktop/combined/TF",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/DEG/dehy_DEG_female.xlsx")


#MALES
DefaultAssay(combined) <- "RNA"
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined)))({
  try({
    ident1 <- paste0(i,"_dehydrated_male")
    ident2 <- paste0(i,"_adlib_male")
    Idents(combined) <- "celltype.groupid.sex"
    deg <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
    deg <-cbind(deg, gene=rownames(deg))
    write.xlsx(deg,file=paste0(output,i,".xlsx"))
  })
})

list_of_files <- list.files(path = "~/Desktop/combined/DEG",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/DEG/dehy_DEG_male.xlsx")


#Ad lib male vs female
DefaultAssay(combined) <- "RNA"
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined)))({
  try({
    ident1 <- paste0(i,"_adlib_male")
    ident2 <- paste0(i,"_adlib_female")
    Idents(combined) <- "celltype.groupid.sex"
    deg <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
    deg <-cbind(deg, gene=rownames(deg))
    write.xlsx(deg,file=paste0(output,i,".xlsx"))
  })
})



Idents(combined) <- "celltype.groupid"
deg <- FindMarkers(combined, ident.1 = "mPTS1_adlib", ident.2="fPTS1_adlib", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/PTS1.xlsx")

deg <- FindMarkers(combined, ident.1 = "mPTS2_adlib", ident.2="fPTS2_adlib", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/PTS2.xlsx")

deg <- FindMarkers(combined, ident.1 = "mPTS3_adlib", ident.2="fPTS3_adlib", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/PTS3.xlsx")

deg <- FindMarkers(combined, ident.1 = "mdTL_adlib", ident.2="fdTL_adlib", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/dTL.xlsx")

setwd("~/Desktop/combined/DEG")
list_of_files <- list.files(path = "~/Desktop/combined/DEG",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/DEG/m_V_female_adlib.xlsx")


#Dehydrated male vs female
DefaultAssay(combined) <- "RNA"
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined)))({
  try({
    ident1 <- paste0(i,"_dehydrated_male")
    ident2 <- paste0(i,"_dehydrated_female")
    Idents(combined) <- "celltype.groupid.sex"
    deg <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
    deg <-cbind(deg, gene=rownames(deg))
    write.xlsx(deg,file=paste0(output,i,".xlsx"))
  })
})


Idents(combined) <- "celltype.groupid"
deg <- FindMarkers(combined, ident.1 = "mPTS1_dehydrated", ident.2="fPTS1_dehydrated", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/PTS1.xlsx")

deg <- FindMarkers(combined, ident.1 = "mPTS2_dehydrated", ident.2="fPTS2_dehydrated", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/PTS2.xlsx")

deg <- FindMarkers(combined, ident.1 = "mPTS3_dehydrated", ident.2="fPTS3_dehydrated", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/PTS3.xlsx")

deg <- FindMarkers(combined, ident.1 = "mdTL_dehydrated", ident.2="fdTL_dehydrated", logfc.threshold = 0, verbose = TRUE, only.pos=FALSE) 
deg <-cbind(deg, gene=rownames(deg))
write.xlsx(deg, file = "~/Desktop/combined/DEG/dTL.xlsx")

setwd("~/Desktop/combined/DEG")
list_of_files <- list.files(path = "~/Desktop/combined/DEG",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/DEG/m_V_female_dehydrated.xlsx")




#DAC for dehydrated vs adlib
DefaultAssay(combined) <- "peaks"
setwd("~/Desktop/combined/TF/")
output <- "~/Desktop/combined/TF/"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated")
    ident2 <- paste0(i,"_adlib")
    Idents(combined) <- "celltype.groupid"
    DAC <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(combined, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
setwd("~/Desktop/combined/TF/")
list_of_files <- list.files(path = "~/Desktop/combined/TF/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <-df %>% dplyr::filter(p_val_adj <0.05)
write.xlsx(new_df, file = "~/Desktop/combined/TF/dehy_DAC.xlsx")


#DAC females
DefaultAssay(combined) <- "peaks"
setwd("~/Desktop/combined/DEG/")
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated_female")
    ident2 <- paste0(i,"_adlib_female")
    Idents(combined) <- "celltype.groupid.sex"
    DAC <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(combined, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
setwd("~/Desktop/combined/DAC/Female")
list_of_files <- list.files(path = "~/Desktop/combined/DAC/Female",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <-df %>% dplyr::filter(p_val_adj <0.05)
write.xlsx(new_df, file = "~/Desktop/combined/DAC/Female/Fem_DAC.xlsx")

#DAC males
DefaultAssay(combined) <- "peaks"
setwd("~/Desktop/combined/DEG/")
output <- "~/Desktop/combined/DEG/"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated_male")
    ident2 <- paste0(i,"_adlib_male")
    Idents(combined) <- "celltype.groupid.sex"
    DAC <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, verbose = TRUE, min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') 
    open <- rownames(DAC)  
    cf <- ClosestFeature(combined, regions = open)
    DAC_final <-(cbind(DAC, gene=cf$gene_name, distance=cf$distance))
    write.xlsx(DAC_final, file=paste0(output,i,".xlsx"))
  })
}
setwd("~/Desktop/combined/DAC/Male")
list_of_files <- list.files(path = "~/Desktop/combined/DAC/Male",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
new_df <-df %>% dplyr::filter(p_val_adj <0.05)
write.xlsx(new_df, file = "~/Desktop/combined/DAC/Male/male_DAC.xlsx")

#################Cell cycle###############################################
#turns out there isn't anything exciting with dehydration and cell cycle.  
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
Idents(combined) <- "celltype.groupid"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#make sentence case to match GEX
s.genes <-str_to_title(s.genes)
g2m.genes <- str_to_title(g2m.genes)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst")
combined <- ScaleData(combined, features = rownames(combined))
combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
combined <- RunPCA(combined, features = c(s.genes, g2m.genes))
DimPlot(combined, split.by = "sex")
RidgePlot(combined, features = "Cdk1", ncol = 2)
RidgePlot(combined, features = c("Cdk1", "Top2a", "Mcm6", "Mki67"), ncol = 2)
VlnPlot(combined, features = "Pcna", ncol = 2)
# view cell cycle scores and phase assignments
head(combined[[]])
#########################################################################################################



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
