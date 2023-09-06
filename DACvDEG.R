library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(openxlsx)
library(readxl)
library(dplyr)

###################T5n insertion##################################################
library(Signac) 
library(Seurat) 
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
library(openxlsx)
set.seed(1234)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(org.Mm.eg.db)

#t5n insertion in Figure S3.
#dehydration separate from ad lib.  Loaded object, formerly called combined, and renamed dehydration.
dehydration <- readRDS("/dehydration.rds")
dehydrated <-subset(dehydration, subset = groupid == "dehydrated")
adlib <-subset(dehydration, subset = groupid == "adlib")
DefaultAssay(dehydrated) <- 'peaks'
DefaultAssay(adlib) <- 'peaks'
Idents(dehydrated) <- "celltype"
Idents(adlib) <- "celltype"

#get dehydrated specific 
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding dac for: ",cluster))
  dehydrated <- seurat_aggregate
  dac <- FindMarkers(dehydrated, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "atac_peak_region_fragments",
                     min.pct = 0.2)
  open <- rownames(dac)  
  cf <- ClosestFeature(dehydrated, regions = open)
  return(cbind(dac, gene=cf$gene_name, distance=cf$distance))}

idents <- levels(dehydrated@meta.data$celltype)
list.cluster.dac <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = dehydrated)})
write.xlsx(list.cluster.dac, file = "/dac.dehydrated.xlsx", sheetName = idents, rowNames = T)

#get adlib specific 
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding dac for: ",cluster))
  adlib <- seurat_aggregate
  dac <- FindMarkers(adllib, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "atac_peak_region_fragments",
                     min.pct = 0.2)
  open <- rownames(dac)  
  cf <- ClosestFeature(adlib, regions = open)
  return(cbind(dac, gene=cf$gene_name, distance=cf$distance))}

idents <- levels(adlib@meta.data$celltype)
list.cluster.dac <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = adlib)})
write.xlsx(list.cluster.dac, file = "/dac.adlib.xlsx", sheetName = idents, rowNames = T)

#open the DAC dehydrated file generated with the DAC_vs_DEG code.
dacfile <- ("/dac.dehydrated.xlsx")
idents <- getSheetNames(dacfile)
list.dac <- lapply(idents, function(x) {
  df <- read.xlsx(dacfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::mutate(celltype = x) 
})

# identify all unique cell-type-specific peaks and filter for log2fc > 0
all_dac <- bind_rows(list.dac) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::select("coord") %>%
  dplyr::distinct()

dac_aver <- AverageExpression(CON, features = all_dac$coord, assays = "peaks")
figxa <- pheatmap::pheatmap(dac_aver[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F,
                            show_rownames = FALSE) #520x630

all_dac.gr <- StringToGRanges(all_dac$coord, sep = c(":","-"))
list.dac.gr <- lapply(seq(list.dac), function(x) {
  df <- list.dac[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dac.gr) <- idents

list.peakAnno <- lapply(list.dac.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)
all.peakAnno <- annotatePeak(all_dac.gr, TxDb = txdb,
                             tssRegion = c(-3000, 3000), verbose = FALSE)

figS1B <- plotAnnoPie(all.peakAnno) #total T5n insertion in the dataset
figS1D <- plotAnnoBar(list.peakAnno) #celltype-specific analysis
figS1F <- plotDistToTSS(list.peakAnno) #Distribution of transcription factor-binding loci relative to TSS

#to do statistics on the distribution export as xlsx
write.xlsx(list.peakAnno, file = ("/data/user/hyndmank/multiomic/dehydration2021/DAC/dehydrated_peakanno.xlsx"))

#same thing but with the adlib files
Idents(adlib) <- "celltype"

dacfile <- ("/dac.adlib.xlsx")
idents <- getSheetNames(dacfile)
list.dac <- lapply(idents, function(x) {
  df <- read.xlsx(dacfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::mutate(celltype = x) 
})

all_dac <- bind_rows(list.dac) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::select("coord") %>%
  dplyr::distinct()

dac_aver <- AverageExpression(KO, features = all_dac$coord, assays = "peaks")
figxb <- pheatmap::pheatmap(dac_aver[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F,
                            show_rownames = FALSE) #520x630

all_dac.gr <- StringToGRanges(all_dac$coord, sep = c(":","-"))
list.dac.gr <- lapply(seq(list.dac), function(x) {
  df <- list.dac[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dac.gr) <- idents
# annotate the list of GRanges dac for each cell type
list.peakAnno <- lapply(list.dac.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)
all.peakAnno <- annotatePeak(all_dac.gr, TxDb = txdb,
                             tssRegion = c(-3000, 3000), verbose = FALSE)

figS1A <- plotAnnoPie(all.peakAnno) 
figS1C <- plotAnnoBar(list.peakAnno) 
figS1E <- plotDistToTSS(list.peakAnno)

#to do statistics on the distribution export as xlsx
write.xlsx(list.peakAnno, file = ("adlib_peakanno.xlsx"))

################DAC versus DEG for the cluster DEG and DACs (cluster of interest versus all other clusters##############
#kept only adjusted p<0.05, import excel where each tab is a different cluster.  make sure tab cluster names are the same in the dac and deg.
#Figure 1B

DAC <- read_excel("~/Desktop/Dehydration/clusterdac.xlsx")
DACdf <- DAC[order(DAC$gene, -abs(DAC$avg_log2FC) ), ]### sort first
DACdf <- DAC %>% distinct(cluster, gene, .keep_all = TRUE) #filter distinct genes in each celltype

DEG <- read_excel("~/Desktop/Dehydration/clusterdeg.xlsx")
DEGdf <- DEG[order(DEG$gene, -abs(DEG$avg_log2FC) ), ]### sort first
DEGdf <- DEG %>% distinct(cluster, gene, .keep_all = TRUE)

#join the two dataframes based upon the Cell (cluster) and gene.
DACDEG <-left_join(DEGdf, DACdf, by=c("cluster","gene"))
write.xlsx(DACDEG, file = ("~/Desktop/Dehydration/DACDEG.xlsx")) 

#downloaded files and removed rows that were not shared between the datasets. Used Prism to do correlations. Plotted in Figure 1B.

#ggplot2 the data if you would like
p <- ggplot(DACDEG, aes(x=avg_log2FC.y, y=avg_log2FC.x, label=gene,xmin=-2.5, xmax = 6.5, ymin =-5, ymax=8))
p + geom_point()

quad_count_CON <- DACDEG %>%
  # Count how many with each combination of X and Y being positive
  dplyr::count(right = avg_log2FC.y > 0, top = avg_log2FC.x > 0) %>%
  # TRUE = 1, FALSE = 0, so these map the TRUE to +1 and FALSE to -1
  mutate(avg_log2FC.y = 2 * (right - 0.5), avg_log2FC.x = 2 * (top - 0.5))
quad_count_CON


