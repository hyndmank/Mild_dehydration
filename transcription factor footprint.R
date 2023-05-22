library(Signac)#1.9.0
library(Seurat)#4.3.0
library(motifmatchr)#1.16.0
library(JASPAR2020)#0.99.10
library(TFBSTools)#1.32.0
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(readxl)
library(openxlsx)
set.seed(1234)
library(chromVAR)
library(Matrix)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(tibble)

DefaultAssay(combined) <- 'peaks'
Idents(combined) <-"celltype"
Idents(combined) <-"celltype.groupid"
Idents(combined) <-"celltype.sex"

# extract position frequency matrices for the motifs:  NCBI tax ids Mus 10090
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# add motif information
combined <- AddMotifs(combined, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
combined

#ChromVAR identifies motifs associated with variability in chromatin accessibility between cells. 
combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(combined) <- 'chromvar'
saveRDS(combined, file = "~/Desktop/Combined/combined.rds")
        
# find motifs with signac.  
DefaultAssay(combined) <- "peaks"
setwd("~/Desktop/combined/TF")
output <- "~/Desktop/combined/TF/"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated")
    ident2 <- paste0(i,"_adlib")
    Idents(combined) <- "celltype.groupid"
    da_peaks <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, only.pos = TRUE,  verbose = TRUE, min.pct = 0.01, test.use = 'LR', latent.vars = 'nCount_peaks')
    top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
    enriched.motifs <- FindMotifs(
      object = combined,
      features = top.da.peak
    )
    write.xlsx(enriched.motifs, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "~/Desktop/combined/TF/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/TF/combined_motif.xlsx")


######################ChromVar###########################################
# create a data frame of motifs and corresponding gene names. Code modified from Muto et al. 2021 PMID: 33850129.
motifs <- combined$peaks@motifs@motif.names
motif_names <- data.frame(genes=unlist(motifs)) %>% rownames_to_column(var="motif") %>% arrange(genes) 

# format the motif_names df so each motif corresponds to a single gene
# the JASPAR database has multiple genes associated with each motif separated by "::"
motif_names <- data.frame(genes=unlist(motifs)) %>% rownames_to_column(var="motif") %>% arrange(genes) 
motif_names <- distinct(motif_names) %>% #subset unique rows
  arrange(genes) %>%
  dplyr::filter(!str_detect(genes, pattern="::"))%>%
  dplyr::filter(!str_detect(genes, pattern="var.")) %>% # remove any genes associated with multiple variants
  as.data.frame()

# fix the one occurrence of EWSR1-FLI1
motif_names$genes <- str_replace(motif_names$genes, pattern="EWSR1-FLI1", replacement="EWSR1")

#Make Motif names lower case to be compatible with gene annotation
motif_names$genes <-paste0(substr(motif_names$genes,1,1),tolower(substr(motif_names$genes,2,nchar(motif_names$genes)))) 
motif.ls = motif_names$motif
gene <- motif_names$genes

#ChromVar Loop.When performing differential testing on the chromVAR z-score, we can set mean.fxn=rowMeans and fc.name="avg_diff" in the FindMarkers() function so that the fold-change calculation computes the average difference in z-score between the groups.  ChromVar doesn't give motif names so I looked them up from peaks motif.

DefaultAssay(combined) <- "chromvar"
setwd("~/Desktop/combined/TF")
output <- "~/Desktop/combined/TF/chrom."
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated")
    ident2 <- paste0(i,"_adlib")
    Idents(combined) <- "celltype.groupid"
    differential.activity <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, features=motif.ls, only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff")
    motifLookup <- rownames(differential.activity)
    motifNames <- sapply(motifLookup, function(x) combined@assays[["peaks"]]@motifs@motif.names [[x]])
    differential.activity <-(cbind(differential.activity, gene = motifNames, motif = motifLookup ))
    write.xlsx(differential.activity, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "~/Desktop/combined/TF/",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/TF/combined_ChromVar.xlsx")

############for motif images###################
DefaultAssay(combined) <- 'peaks'
da_peaks <- FindMarkers(
  object = combined,
  ident.1 = 'CDPC_dehydrated',
  ident.2 = 'CDPC_adlib',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks. 
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
enriched.motifs <- FindMotifs(
  object = combined,
  features = top.da.peak
)

#Crem Motif 400X193 Figure 5G.
MotifPlot(
  object = combined,
  motifs = "CREM"
)

###################################sex differences#################
DefaultAssay(combined) <- "peaks"
setwd("~/Desktop/combined/TF/")
output <- "~/Desktop/combined/TF/female."
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated_female")
    ident2 <- paste0(i,"_adlib_female")
    Idents(combined) <- "celltype.groupid.sex"
    da_peaks <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, only.pos = TRUE,  verbose = TRUE, min.pct = 0.01, test.use = 'LR', latent.vars = 'nCount_peaks')
    top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
    enriched.motifs <- FindMotifs(
      object = combined,
      features = top.da.peak
    )
    write.xlsx(enriched.motifs, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "~/Desktop/combined/TF",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/TF/female_motif.xlsx")

#MALES
DefaultAssay(combined) <- "peaks"
setwd("~/Desktop/combined/TF")
output <- "~/Desktop/combined/TF/male"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated_male")
    ident2 <- paste0(i,"_adlib_male")
    Idents(combined) <- "celltype.groupid.sex"
    da_peaks <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, only.pos = TRUE, verbose = TRUE, min.pct = 0.01, test.use = 'LR', latent.vars = 'nCount_peaks')
    top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
    enriched.motifs <- FindMotifs(
      object = combined,
      features = top.da.peak
    )
    write.xlsx(enriched.motifs, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "~/Desktop/combined/TF",
                            full.names = FALSE)
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/TF/male_motif.xlsx")


####################ChromVar Loop###################################
DefaultAssay(combined) <- "chromvar"
setwd("~/Desktop/combined/TF")
output <- "~/Desktop/combined/TF/"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated_female")
    ident2 <- paste0(i,"_adlib_female")
    Idents(combined) <- "celltype.groupid.sex"
    differential.activity <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, features=motif.ls, only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff")
    motifLookup <- rownames(differential.activity)
    motifNames <- sapply(motifLookup, function(x) combined@assays[["peaks"]]@motifs@motif.names [[x]])
    differential.activity <-(cbind(differential.activity, gene = motifNames, motif = motifLookup ))
    write.xlsx(differential.activity, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "~/Desktop/combined/TF/",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/TF/female_ChromVar.xlsx")

#Males
DefaultAssay(combined) <- "chromvar"
setwd("~/Desktop/combined/TF")
output <- "~/Desktop/combined/TF/"
Idents(combined) <- "celltype"
for (i in (levels(combined))){
  try({
    ident1 <- paste0(i,"_dehydrated_male")
    ident2 <- paste0(i,"_adlib_male")
    Idents(combined) <- "celltype.groupid.sex"
    differential.activity <- FindMarkers(combined, ident.1 = ident1, ident.2=ident2, features=motif.ls, only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff")
    motifLookup <- rownames(differential.activity)
    motifNames <- sapply(motifLookup, function(x) combined@assays[["peaks"]]@motifs@motif.names [[x]])
    differential.activity <-(cbind(differential.activity, gene = motifNames, motif = motifLookup ))
    write.xlsx(differential.activity, file=paste0(output,i,".xlsx"))
  })
}

list_of_files <- list.files(path = "~/Desktop/combined/TF/",
                            full.names = FALSE)
list_of_files
df <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read.xlsx(.x,  colNames = TRUE), .id = "sheet_name")  
write.xlsx(df, file = "~/Desktop/combined/TF/male_ChromVar.xlsx")

#export Z-score matrtix this is a large file
zscore.matrix <-as.matrix(combined [["chromvar"]]@data)
write.csv(zscore.matrix, file = "~/Desktop/combined/TF/zscore.csv")


#chromVar images. Need to use Motif name. MA0609.2 is Crem.
DefaultAssay(combined) <- "chromvar"
p2 <- FeaturePlot(
  object = combined,
  features = "MA0609.2",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = 'umap.rna',
  
)
p2

p3 <- VlnPlot(
  object = combined,
  features = "MA0609.2",
  pt.size = 0.1,
  assay = 'peaks', split.by = "groupid", idents = 'CDPC'
 
)
p3

#########################TF foot printing###############################
# gather the footprinting information for sets of motifs. 
DefaultAssay(combined) <- 'peaks'
Idents(combined) <-"celltype"
CDPC<-subset(x = combined, idents = "CDPC")
DefaultAssay(combined) <- "peaks"
CDPC <- Footprint(
  object = CDPC,
  motif.name = c("CREB1", "CREM", "NFAT5"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)
Figure5H <- PlotFootprint(CDPC, features = c("CREM"), group.by = "groupid") 
Figure5H + patchwork::plot_layout(ncol = 1) + NoLegend()#624x447

#######################correlating DEG with Motifs#########################
#ran the DEGs and Motifs separate.  Modified code from Gerhardt et al. 2023 PMID: 36735940
#Read in list with differentially expressed genes
DEG <- read.xlsx("~/Desktop/combined/dehy_DEG.xlsx")

#Consider only significant genes
DEG <- DEG[which(DEG$p_val_adj<=0.05),]
colnames(DEG) <- c(paste0(rep("GEX." ,7),colnames(DEG)[1:7]), "cluster", "gene")
DEG$gene <- toupper(DEG$gene)

#Read in list with differentially active motifs
Motif <- read.xlsx("~/Desktop/combined/ChromVar.xlsx")
Motif <- Motif[which(Motif$p_val_adj<=0.05),]
colnames(Motif) <- c(paste0(rep("motif." ,8),colnames(Motif)[1:8]), "cluster","gene")
Motif$gene <- toupper(Motif$gene)

#Retain only differentially active motifs that are also differentially expressed
OVERLAP <- data.frame("cluster"=c(),"gene"=c(),"motif"=c(),"motif.avg_diff"=c(), "motif.p_val"=c(),"motif.p_val_adj"=c(),"motif.pct.1"=c(),"motif.pct.2"=c(), 
                      "GEX.avg_log2FC"=c(),"GEX.p_val"=c(),"GEX.p_val_adj"=c(), "GEX.pct.1"=c(),"GEX.pct.2"=c())

Celltypes <- unique(Motif$cluster)

for(i in 1:length(Celltypes)){
  DEG_celltype <- DEG[which(DEG$cluster == Celltypes[[i]]),]
  Motif_celltype <- Motif[which(Motif$cluster== Celltypes[[i]]),]
  Motif_overlap <- Motif_celltype[which(Motif_celltype$gene%in%DEG_celltype$gene),]
  DEG_overlap <- DEG_celltype[which(DEG_celltype$gene%in%Motif_overlap$gene),]
  myDF <- merge(Motif_overlap,DEG_overlap,by=c("gene","cluster"))
  OVERLAP <- rbind(OVERLAP,myDF)
}

#Calculate a Motif_GEX_score
OVERLAP$Motif_GEX_score <- OVERLAP$motif.avg_diff* OVERLAP$GEX.avg_log2FC
write.xlsx(OVERLAP, file = "~/Desktop/combined/TF/motif_chromvar.xlsx")

#Heatmap for the TF activity and Expression
#took the zscores and log2FC and made heatmap in prism. Figure 5F.


