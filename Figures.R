library(Seurat)#v4.3
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

# first compute the GC content for each peak
combined <- RegionStats(combined, genome = BSgenome.Mmusculus.UCSC.mm10)
saveRDS(combined, file = "~/Desktop/Dehydration/combined.rds")

# link peaks to genes of interest
combined <- LinkPeaks(
  object = combined,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Aqp2", "Aqp3", "Aqp4", "Aqp1", "Aqp5", "Avpr2", "Pde10a", "Atp1a1", "Atp1b1", "Grem2", "Btc", "Crem", "Grem1", "Nfat5")
)

#Visualizations Used in the publication.
DimPlot(object = combined, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = combined, label = FALSE, repel = TRUE) + NoLegend() #1030x831 Figure 1A.
DimPlot(object = combined, label = FALSE, repel = TRUE, split.by = "groupid") + NoLegend() #1204x831 Figure 1A insert.
DimPlot(object = combined, label = FALSE, repel = TRUE, group.by = "sex") + NoLegend() #1204x831 Figure S3.

#AQP2 700x400
FeaturePlot(combined, features = "Aqp2", split.by = "groupid", min.cutoff = 0, max.cutoff =3) #Figure 2A
FeaturePlot(combined, features = "Aqp2", min.cutoff = 0, max.cutoff =3)
#RNA summary 250x400 Figure 2B.
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#chromatin summary 250x400 Figure 2C
DefaultAssay(combined) <- "Gene"
VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "CDPC",
  slot = "data",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Figure2D
DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp2",
    idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000,) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#AQP2 CNT Figure 2E
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "CNT1",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#chromatin summary 250x400
VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "CNT1",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 
  
#IMCD Figure 2F
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"

VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "IMCD",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "IMCD",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Figure 2G is normalized counts from the BulkRNA-seq and plotted in prism.

#AQP3 250x400 figure 3A, B
FeaturePlot(combined, features = "Aqp3", split.by = "groupid", min.cutoff = 0, max.cutoff =3)
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp3",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#chromatin summary 250x400 Figure 3C
VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp3",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Figure 3D
DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp3",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 6000,
  extend.downstream = 5000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#AQP4 Figure 3F
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp4",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp4",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Grem2 Figure 4A
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Grem2", split.by = "groupid", min.cutoff = 0, max.cutoff =3)

#Figure 4B, C
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Grem2",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Grem2",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Figure 4D
DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Grem2",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#Pde10a Figure 5.
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Pde10a", split.by = "groupid", min.cutoff = 0, max.cutoff =3)

VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Pde10a",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Pde10a",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Pde10a",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 20000,
  extend.downstream = 20000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#DCT 500x400  Figure 5E, F
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Pde10a",
  pt.size = 0.5,
  idents = c("DCT1", "DCT2"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Pde10a",
  pt.size = 0.5,
  idents = c("DCT1", "DCT2"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Podocyte 250x400 Figure 5G
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Pde10a",
  pt.size = 0.5,
  idents = "Podocyte",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Pde10a",
  pt.size = 0.5,
  idents = "Podocyte",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#CREM Figure 6A-D
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Crem", split.by = "groupid", min.cutoff = 0, max.cutoff =3)

VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Crem",
  pt.size = 0.5,
  idents = 'CDPC',
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue"))+ NoLegend()  & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Crem",
  pt.size = 0.5,
  idents = 'CDPC',
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue"))+ NoLegend()  & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Crem",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000, annotation = TRUE, links = TRUE) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#NFAT5. 500x250 Figure S2.
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Nfat5", split.by = "groupid", min.cutoff = 0, max.cutoff =3)

VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Nfat5",
  pt.size = 0.5,
  idents = c("CNT1", "CDPC", "IMCD", "VR1"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Aqp1 Figure S4
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Aqp1", split.by = "groupid", min.cutoff = 0, max.cutoff =3)

#700x400
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp1",
  pt.size = 0.5,
  idents = c("mPTS1", "mPTS2", "mPTS3", "mdTL", "dTL", "VR1", "VR2", "VR3"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue"))  & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp1",
  pt.size = 0.5,
  idents = c("mPTS1", "mPTS2", "mPTS3", "mdTL", "dTL", "VR1", "VR2", "VR3"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#500 X 400
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp1",
  pt.size = 0.5,
  idents = c("fPTS1", "fPTS2", "fPTS3", "fdTL"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue"))+ NoLegend()  & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp1",
  pt.size = 0.5,
  idents = c("fPTS1", "fPTS2", "fPTS3", "fdTL"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 
