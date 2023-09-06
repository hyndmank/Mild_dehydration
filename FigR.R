################ Peak-gene association testing with FIGR ##############
#R version 4.2.2
library(dplyr)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)
library(writexl)
set.seed(123)

######################### All clusters #########################
### Load Seurat object
load("dehydration.RData")

############################# Ad lib ###########################
#Subset Adlib clusters from Seurat object
Adlib <- subset(dehydration, subset = groupid == "adlib")
save(Adlib,file="Adlib.RData")

#Create RangedSummarizedExperiment object for ATAC from Adlib Seurat object
Adlib_ATAC.se <- SummarizedExperiment(assays=Adlib[["peaks"]]@counts,  
                                      rowRanges=Adlib[["peaks"]]@ranges)
assayNames(Adlib_ATAC.se) <- "counts"

save(Adlib_ATAC.se,file="Adlib_ATAC_se.RData")

#Get scRNA-seq matrix from Seurat object
Adlib_RNAmat <- Adlib[["RNA"]]@counts

dim(Adlib_ATAC.se) # Peaks x Cells
dim(Adlib_RNAmat) # Genes x Cells

#Remove genes with zero expression across all cells
Adlib_RNAmat <- Adlib_RNAmat[Matrix::rowSums(Adlib_RNAmat)!=0,]
dim(Adlib_RNAmat) # Genes x Cells
save(Adlib_RNAmat,file="Adlib_RNAmat.RData")

#Peak-gene association testing
Adlib_cisCorr <- FigR::runGenePeakcorr(ATAC.se = Adlib_ATAC.se,
                                       RNAmat = Adlib_RNAmat,
                                       genome = "mm10", 
                                       nCores = 10,
                                       keepPosCorOnly = TRUE,
                                       p.cut = NULL, # Set this to NULL and we can filter later
                                       n_bg = 100)
head(Adlib_cisCorr)
save(Adlib_cisCorr,file="~Desktop/All_Adlib_cisCorr.RData")
write_xlsx(Adlib_cisCorr,"~Desktop/All_Adlib_cisCorr.xlsx")

#Filter peak-gene pair with significant association
Adlib_cisCorr.filt <- Adlib_cisCorr %>% filter(pvalZ <= 0.05)
write_xlsx(Adlib_cisCorr.filt,"~Desktop/All_Adlib_cisCorr_filtered.xlsx")


############################## Dehydrated #############################################
#Subset Dehydrated clusters from Seurat object
Deh <- subset(dehydration, subset = groupid == "dehydrated")
save(Deh,file="Deh.RData")

#Create RangedSummarizedExperiment object for ATAC from Dehydrated Seurat object

Deh_ATAC.se <- SummarizedExperiment(assays=Deh[["peaks"]]@counts,  
                                    rowRanges=Deh[["peaks"]]@ranges)
assayNames(Deh_ATAC.se) <- "counts"
save(Deh_ATAC.se,file="~Desktop/Deh_ATAC_se.RData")

#Get scRNA-seq matrix from Seurat object
Deh_RNAmat <- Deh[["RNA"]]@counts
dim(Deh_ATAC.se) # Peaks x Cells
dim(Deh_RNAmat) # Genes x Cells

#Remove genes with zero expression across all cells
Deh_RNAmat <- Deh_RNAmat[Matrix::rowSums(Deh_RNAmat)!=0,]
save(Deh_RNAmat,file="~Desktop/Deh_RNAmat.RData")
dim(Deh_RNAmat) # Genes x Cells

#Peak-gene association testing
Deh_cisCorr <- FigR::runGenePeakcorr(ATAC.se = Deh_ATAC.se,
                                     RNAmat = Deh_RNAmat,
                                     genome = "mm10",  
                                     nCores = 10,
                                     keepPosCorOnly = TRUE,
                                     p.cut = NULL, # Set this to NULL and we can filter later
                                     n_bg = 100)
head(Deh_cisCorr)
save(Deh_cisCorr,file="~Desktop/All_Deh_cisCorr.RData")
write_xlsx(Deh_cisCorr,"~Desktop/All_Deh_cisCorr.xlsx")

#Filter peak-gene pairs with significant association
Deh_cisCorr.filt <- Deh_cisCorr %>% filter(pvalZ <= 0.05)
write_xlsx(Deh_cisCorr.filt,"~Desktop/All_Deh_cisCorr_filtered.xlsx")


########################## Principle cell ############################
#Subset PC Seurat
PC <- subset(dehydration, idents = c("CDPC", "CNT-PC"))
save(PC,file="~Desktop/PC.RData")

#Subset PC ad lib 
PC_Adlib <- subset(PC, subset = groupid == "adlib")
save(PC_Adlib,file="~Desktop/PC_Adlib.RData")

#Subset PC dehydrated
PC_Deh <- subset(PC, subset = groupid == "dehydrated")
save(PC_Deh,file="~Desktop/PC_Deh.RData")

######################### PC Ad lib ##################################

#Create RangedSummarizedExperiment object for ATAC from PC_Adlib Seurat object
PC_Adlib_ATAC.se <- SummarizedExperiment(assays=PC_Adlib[["peaks"]]@counts,  
                                         rowRanges=PC_Adlib[["peaks"]]@ranges)
assayNames(PC_Adlib_ATAC.se) <- "counts"
save(PC_Adlib_ATAC.se,file="PC_Adlib_ATAC_se.RData")

#Get scRNA-seq matrix from PC_Adlib Seurat object
PC_Adlib_RNAmat <- PC_Adlib[["RNA"]]@counts
dim(PC_Adlib_ATAC.se) # Peaks x Cells
dim(PC_Adlib_RNAmat) # Genes x Cells

#Remove genes with zero expression across all cells
PC_Adlib_RNAmat <- PC_Adlib_RNAmat[Matrix::rowSums(PC_Adlib_RNAmat)!=0,]
save(PC_Adlib_RNAmat,file="PC_Adlib_RNAmat.RData")

#Peak-gene association testing
PC_Adlib_cisCorr <- FigR::runGenePeakcorr(ATAC.se = PC_Adlib_ATAC.se,
                                          RNAmat = PC_Adlib_RNAmat,
                                          genome = "mm10", # One of hg19, mm10 or hg38 
                                          keepPosCorOnly = TRUE,
                                          nCores = 10,
                                          p.cut = NULL, # Set this to NULL and we can filter later
                                          n_bg = 100)
head(PC_Adlib_cisCorr)
save(PC_Adlib_cisCorr,file="~Desktop/PC_Adlib_cisCorr.RData")
write_xlsx(PC_Adlib_cisCorr,"~Desktop/PC_Adlib_cisCorr.xlsx")

#Filter peak-gene pairs with significant association
PC_Adlib_cisCorr.filt <- PC_Adlib_cisCorr %>% filter(pvalZ <= 0.05)
save(PC_Adlib_cisCorr.filt, file="~Desktop/PC_Adlib_cisCorr_filt.RData")
write_xlsx(PC_Adlib_cisCorr.filt,"~Desktop/PC_Adlib_cisCorr.filtered.xlsx")

############################## PC Dehydrated ########################################
#### Build RangedSummarizedExperiment object for ATAC from PC_Deh Seurat object
PC_Deh_ATAC.se <- SummarizedExperiment(assays=PC_Deh[["peaks"]]@counts,  
                                       rowRanges=PC_Deh[["peaks"]]@ranges)
assayNames(PC_Deh_ATAC.se) <- "counts"
save(PC_Deh_ATAC.se,file="PC_Deh_ATAC_se.RData")

#Get scRNA-seq matrix from PC_Deh Seurat object
PC_Deh_RNAmat <- PC_Deh[["RNA"]]@counts

dim(PC_Deh_ATAC.se) # Peaks x Cells
dim(PC_Deh_RNAmat) # Genes x Cells

#Remove genes with zero expression across all cells
PC_Deh_RNAmat <- PC_Deh_RNAmat[Matrix::rowSums(PC_Deh_RNAmat)!=0,]
save(PC_Deh_RNAmat,file="PC_Deh_RNAmat.RData")

#Peak-gene association testing
PC_Deh_cisCorr <- FigR::runGenePeakcorr(ATAC.se = PC_Deh_ATAC.se,
                                        RNAmat = PC_Deh_RNAmat,
                                        genome = "mm10", # One of hg19, mm10 or hg38 
                                        keepPosCorOnly = TRUE,
                                        nCores = 10,
                                        p.cut = NULL, # Set this to NULL and we can filter later
                                        n_bg = 100)
head(PC_Deh_cisCorr)
save(PC_Deh_cisCorr,file="~Desktop/PC_Deh_cisCorr.RData")
write_xlsx(PC_Deh_cisCorr,"~Desktop/PC_Deh_cisCorr.xlsx")

#Filter peak-gene pairs with significant association
PC_Deh_cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
save(PC_Deh_cisCorr.filt, file="~Desktop/PC_Deh_cisCorr_filt.RData")
write_xlsx(PC_Deh_cisCorr.filt,"~Desktop/PC_Deh_cisCorr.filtered.xlsx")

############################# Proximal tubules ################################
#Subset PT Seurat
PT <- subset(dehydration, idents = c("fPTS1", "fPTS2", "fPTS2", "mPTS1", "mPTS2", "mPTS3"))
save(PT,file="~Desktop/PT.RData")

#Subset PT ad lib 
PT_Adlib <- subset(PT, subset = groupid == "adlib")
save(PT_Adlib,file="~Desktop/PT_Adlib.RData")

#Subset PT dehydrated
PT_Deh <- subset(PT, subset = groupid == "dehydrated")
save(PT_Deh,file="~Desktop/PT_Deh.RData")

######################## PT Ad lib ################################
#Create RangedSummarizedExperiment object for ATAC from PT_Adlib Seurat object
PT_Adlib_ATAC.se <- SummarizedExperiment(assays=PT_Adlib[["peaks"]]@counts,  
                                         rowRanges=PT_Adlib[["peaks"]]@ranges)
assayNames(PT_Adlib_ATAC.se) <- "counts"
save(PT_Adlib_ATAC.se,file="PT_Adlib_ATAC_se.RData")

#Get scRNA-seq matrix from PT_Adlib Seurat object
PT_Adlib_RNAmat <- PT_Adlib[["RNA"]]@counts
dim(PT_Adlib_ATAC.se) # Peaks x Cells
dim(PT_Adlib_RNAmat) # Genes x Cells

#Remove genes with zero expression across all cells
PT_Adlib_RNAmat <- PT_Adlib_RNAmat[Matrix::rowSums(PT_Adlib_RNAmat)!=0,]
save(PT_Adlib_RNAmat,file="PT_Adlib_RNAmat.RData")

#Peak-gene association testing
PT_Adlib_cisCorr <- FigR::runGenePeakcorr(ATAC.se = PT_Adlib_ATAC.se,
                                          RNAmat = PT_Adlib_RNAmat,
                                          genome = "mm10", # One of hg19, mm10 or hg38 
                                          keepPosCorOnly = TRUE,
                                          nCores = 10,
                                          p.cut = NULL, # Set this to NULL and we can filter later
                                          n_bg = 100)
head(PT_Adlib_cisCorr)
save(PT_Adlib_cisCorr,file="~Desktop/PT_Adlib_cisCorr.RData")
write_xlsx(PT_Adlib_cisCorr,"~Desktop/PT_Adlib_cisCorr.xlsx")

#Filter peak-gene pairs with significant association
PT_Adlib_cisCorr.filt <- PT_Adlib_cisCorr %>% filter(pvalZ <= 0.05)
save(PT_Adlib_cisCorr.filt, file="~Desktop/PT_Adlib_cisCorr_filt.RData")
write_xlsx(PT_Adlib_cisCorr.filt,"~Desktop/PT_Adlib_cisCorr.filtered.xlsx")

###################### PT Dehydrated ########################
#Build RangedSummarizedExperiment object for ATAC from PT_Deh Seurat object
PT_Deh_ATAC.se <- SummarizedExperiment(assays=PT_Deh[["peaks"]]@counts,  
                                       rowRanges=PT_Deh[["peaks"]]@ranges)
assayNames(PT_Deh_ATAC.se) <- "counts"
save(PT_Deh_ATAC.se,file="PT_Deh_ATAC_se.RData")

#Get scRNA-seq matrix from PT_Deh Seurat object
PT_Deh_RNAmat <- PT_Deh[["RNA"]]@counts
dim(PT_Deh_ATAC.se) # Peaks x Cells
dim(PT_Deh_RNAmat) # Genes x Cells

#Remove genes with zero expression across all cells
PT_Deh_RNAmat <- PT_Deh_RNAmat[Matrix::rowSums(PT_Deh_RNAmat)!=0,]
save(PT_Deh_RNAmat,file="PT_Deh_RNAmat.RData")

#Peak-gene association testing
PT_Deh_cisCorr <- FigR::runGenePeakcorr(ATAC.se = PT_Deh_ATAC.se,
                                        RNAmat = PT_Deh_RNAmat,
                                        genome = "mm10", # One of hg19, mm10 or hg38 
                                        keepPosCorOnly = TRUE,
                                        nCores = 10,
                                        p.cut = NULL, # Set this to NULL and we can filter later
                                        n_bg = 100)

head(PT_Deh_cisCorr)
save(PT_Deh_cisCorr,file="E:/Van/R/PT_Deh_cisCorr.RData")
write_xlsx(PT_Deh_cisCorr,"E:\\Van\\R\\PT_Deh_cisCorr.xlsx")

#Filter peak-gene pairs with significant association
PT_Deh_cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
save(PT_Deh_cisCorr.filt, file="E:/Van/R/PT_Deh_cisCorr_filt.RData")
write_xlsx(PT_Deh_cisCorr.filt,"E:\\Van\\R\\PT_Deh_cisCorr.filtered.xlsx")


################### Heatmaps of top genes with significant peak-gene association ################################

#Filtered gene-peak pairs tables to get protein-coding genes only (using EnsDb.Mmusculus.v79 gene annotation) 
#Filter only p< 0.05
#Choose top 50 genes with the highest correlation coefficient (rObs)

#Draw heatmaps comparing Ad lib vs. Dehydrated 
library(ggplot2)
library("scales")

#All clusters. Figure 1C.
mat <- read.csv("All_AdlibvsDeh_Top50.csv")
mycol <- c("yellow","red","red4")
ggplot(data = mat, aes(x = Treatment, y = Gene, fill = rObs)) +
  geom_tile() +
  scale_fill_gradientn(colours = mycol) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Principal cells. Figure S4A
mat <- read.csv("PC_AdlibvsDeh_Top50.csv")
ggplot(data = mat, aes(x = Treatment, y = Gene, fill = rObs)) +
  geom_tile() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_gradient2(low="yellow", mid="orange", high="red", #colors in the scale
                       midpoint=0.4, #same midpoint for plots 
                       na.value = "grey50",   
                       breaks=seq(0.2,1,0.2),#breaks in the scale bar
                       limits=c(0.2, 1),
                       guide = "colourbar") 

#Proximal tubules. Figure S4B
mat <- read.csv("PT_AdlibvsDeh_Top50.csv")
ggplot(data = mat, aes(x = Treatment, y = Gene, fill = rObs)) +
  geom_tile() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_gradient2(low="yellow", mid="orange", high="red", #colors in the scale
                       midpoint=0.4, #same midpoint for plots 
                       na.value = "grey50",   
                       breaks=seq(0.2,1,0.2),#breaks in the scale bar
                       limits=c(0.2, 1),
                       guide = "colourbar") 
