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

# link peaks to genes
combined <- LinkPeaks(
  object = combined,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Aqp2", "Aqp3", "Aqp4", "Aqp1", "Aqp5", "Avpr2", "Pde10a", "Atp1a1", "Atp1b1", "Grem2", "Btc", "Crem", "Grem1", "Nfat5")
)

#Visualizations Used in the publication.
DimPlot(object = combined, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = combined, label = FALSE, repel = TRUE) + NoLegend() #1030x831
DimPlot(object = combined, label = TRUE, repel = TRUE) 
DimPlot(object = combined, label = FALSE, repel = TRUE) + NoLegend() #1205x831
DimPlot(object = combined, label = FALSE, repel = TRUE, group.by = "groupid") #1204x831
DimPlot(object = combined, label = FALSE, repel = TRUE, split.by = "groupid") + NoLegend() #1204x831
DimPlot(object = combined, label = FALSE, repel = TRUE, split.by = "sex") + NoLegend() #1204x831
DimPlot(object = combined, label = FALSE, repel = TRUE, group.by = "sex") + NoLegend() #1204x831

#AQP2 700x400
FeaturePlot(combined, features = "Aqp2", split.by = "groupid", min.cutoff = 0, max.cutoff =3) #
FeaturePlot(combined, features = "Aqp2", min.cutoff = 0, max.cutoff =3)
#RNA summary 250x400
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#chromatin summary 250x400
DefaultAssay(combined) <- "Gene"
VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp2",
  pt.size = 0.5,
  idents = "CDPC",
  slot = "data",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp2",
    idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000,) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#AQP2 CNT
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

#coverageplot 863x300
DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp2",
  idents = c("CNT2_adlib", "CNT2_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000)& scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))
  
#IMCD
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


#AQP3 250x400
FeaturePlot(combined, features = "Aqp3", split.by = "groupid", min.cutoff = 0, max.cutoff =3)
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp3",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#chromatin summary 250x400
VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp3",
  pt.size = 0.5,
  idents = "CDPC",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp3",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 6000,
  extend.downstream = 5000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#IMCD aqp3
VlnPlot(
object = combined,
assay = 'RNA',
features = "Aqp3",
pt.size = 0.5,
idents = "CNT2",
split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 


#AQP4
FeaturePlot(combined, features = "Aqp4", split.by = "groupid", min.cutoff = 0, max.cutoff =3)
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

#AQP4
VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Aqp4",
  pt.size = 0.5,
  idents = "IMCD",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Aqp4",
  pt.size = 0.5,
  idents = "IMCD",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

#Grem2
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Grem2", split.by = "groupid", min.cutoff = 0, max.cutoff =3)
FeaturePlot(combined, features = "Grem2", split.by = "sex", min.cutoff = 0, max.cutoff =3)
FeaturePlot(combined, features = "Grem2", min.cutoff = 0, max.cutoff =3)

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

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Grem2",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#Pde10a
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

#DCT 500x400 
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

#Podocyte 250x400
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

#NFAT5. this was downregulated in a few cell types.500x250
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "celltype"
FeaturePlot(combined, features = "Nfat5", split.by = "groupid", min.cutoff = 0, max.cutoff =3)
FeaturePlot(combined, features = "Nfat5", min.cutoff = 0, max.cutoff =3)

VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Nfat5",
  pt.size = 0.5,
  idents = c("CNT1", "CDPC", "IMCD", "VR1"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

VlnPlot(
  object = combined,
  assay = 'Gene',
  features = "Nfat5",
  pt.size = 0.5,
  idents = "fPTS1",
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Nfat5",
  idents = c("IMCD_adlib", "IMCD_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 20000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#Aqp1
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

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp1",
  idents = c("mPTS3_adlib", "mPTS3_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))

#CREM
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

VlnPlot(
  object = combined,
  assay = 'RNA',
  features = "Crem",
  pt.size = 0.5,
  idents = c("mPTS1", "mPTS2", "mPTS3", "mdTL", "dTL", "VR1", "VR2", "VR3"),
  split.by = "groupid",cols=c("darkolivegreen3", "deepskyblue")) + NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_blank(), axis.title.y = element_blank()) 

DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Crem",
  idents = c("CDPC_adlib", "CDPC_dehydrated"),
  extend.upstream = 10000,
  extend.downstream = 10000, annotation = TRUE, links = TRUE) & scale_fill_manual(values = c("darkolivegreen3", "deepskyblue"))


#AQP5
DefaultAssay(combined) <- "peaks"
Idents(combined) <- "celltype.groupid"
CoveragePlot(
  object = combined,
  region = "Aqp5", idents = c("CDPC_adlib", "CDPC_dehydrated"),
    extend.upstream = 10000,
  extend.downstream = 10000) 

idents.plot <- c("CDPC_adlib", "CDPC_dehydrated")
idents.plot.m <- c("mPTS1_adlib", "mPTS1_dehydrated", "mPTS2_adlib", "mPTS2_dehydrated", "mPTS3_adlib", "mPTS3_dehydrated")
idents.plot.f <- c("fPTS1_adlib", "fPTS1_dehydrated", "fPTS2_adlib", "fPTS2_dehydrated", "fPTS3_adlib", "fPTS3_dehydrated")


p1 <- CoveragePlot(
  object = combined,
  region = "Aqp2",
  features = "Aqp2",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000
)
p1
p2 <- CoveragePlot(
  object = combined,
  region = "Aqp3",
  features = "Aqp3",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)

p3 <- CoveragePlot(
  object = combined,
  region = "Aqp1",
  features = "Aqp1",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000
)

p4 <- CoveragePlot(
  object = combined,
  region = "Aqp4",
  features = "Aqp4",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)
p5 <- CoveragePlot(
  object = combined,
  region = "Avpr2",
  features = "Avpr2",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000, annotation = TRUE)


p5 + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

patchwork::wrap_plots(p3, p4, ncol = 1)
patchwork::wrap_plots(p1,p2, p4, ncol = 1)
p5

p6 <- CoveragePlot(
  object = combined,
  region = "Pde10a",
  features = "Pde10a",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000, annotation = TRUE)
p6 #800x625
p7 <- CoveragePlot(
  object = combined,
  region = "Atp1a1",
  features = "Atp1a1",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000, annotation = TRUE)
p7

p8 <- CoveragePlot(
  object = combined,
  region = "Atp1b1",
  features = "Atp1b1",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000, annotation = TRUE)
p8

p9 <- CoveragePlot(
  object = combined,
  region = "Aqp1",
  features = "Aqp1",
  expression.assay = "RNA",
  idents = idents.plot.m,
  extend.upstream = 8000,
  extend.downstream = 5000, annotation = TRUE)
p9


p10 <- CoveragePlot(
  object = combined,
  region = "Aqp1",
  features = "Aqp1",
  expression.assay = "RNA",
  idents = idents.plot.f,
  extend.upstream = 8000,
  extend.downstream = 5000, annotation = TRUE)
p10
