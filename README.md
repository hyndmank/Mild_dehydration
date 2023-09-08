
# <p align="center"><b>Mild dehydration effects on the murine kidney single nucleus transcriptome and chromatin accessibility.</hx></b></br>
<p align="center"> Nha Van Huynh and Kelly A. Hyndman  <br>
University of Alabama at Birmingham <br></c>


10X Genomics single nucleus multiome-Seq and bulk kidney RNA-seq.  Male and female C57bl/6J mice (10-11 weeks old) were either left with water  (ad lib) or water deprived (dehydrated) for 9-11 hours.  Euthanasia and sample collection occurred +/-1 hour of midnight.

Half kidneys were used for either nuclei isolation and snRNA-seq and snATAC-seq, or bulk RNA-seq (validation cohort).
For the sn-multiome:
630	Female, Ad lib: 636	Female, water deprived: 638	Male, water deprived: 641	Male, Ad lib.  Raw and analyzed files can be found at Geo GSE232662.<br>

CellRanger ARC2.0.0 was used for initial counts and mapping,and then [Seurat](https://satijalab.org/seurat/) and [Signac](https://stuartlab.org/signac/) were used for final analyses. Our code was modified taking into consideration the code by the [Humphreys lab](https://github.com/p4rkerw/Muto_Wilson_NComm_2020/tree/master) and the [McMahon lab](https://github.com/lmsgerhardt/Gerhardt_JASN_2022).  We also used the [FIGR package](https://github.com/buenrostrolab/FigR).<br>

For the bulk RNA-seq cohort:  FA1	female ad lilb 1, FA2	female ad lilb 2, FA3	female ad lilb 3, FA4	female ad lilb 4, FA5	female ad lilb 5, FA6	female ad lilb 6, FD1	female dehydrated 1, FD2	female dehydrated 2, FD4	female dehydrated 3, FD5	female dehydrated 4, FD6	female dehydrated 5, MA1	male ad lib 1, MA2	male ad lib 2, MA3	male ad lib 3, MA4	male ad lib 4, MA5	male ad lib 5, MA6	male ad lib 6, MD1	male dehydrated 1, MD2	male dehydrated 2, MD3	male dehydrated 3, MD4	male dehydrated 4, MD5	male dehydrated 5, MD6	male dehydrated 6, MD7	male dehydrated 7.  Raw and analyzed files can be found at Geo GSE232792.<br>

Fastq files were trimmed with [Trim_Galore](https://github.com/FelixKrueger/TrimGalore), mapped with [Star](https://github.com/alexdobin/STAR), and counts with featurecounts function of [subread](https://subread.sourceforge.net).  All normalization and differentially expressed genes were calucated using the [Deseq2 package](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).


Publication: Mild dehydration effects on the murine kidney single nucleus transcriptome and chromatin accessibility.  Nha Van Huynh, Cassidy Rehage, and Kelly A. Hyndman, University of Alabama at Birmingham, Department of Medicine, Division of Nephrology, Section of Cardio-Renal Physiology and Medicine, Birmingham, AL 35233.  Under Review with American Journal of Physiology: Renal Physiology.
