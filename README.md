# GapCNV
GapCNV: A Comprehensive Framework for CNV Detection in Low-input DNA Sequencing Data

# Author
Xuanxuan Yu, Fei Qin, Shiwei Liu, Noah Brown, Xiao Fan, Qing Lu, Guoshuai Cai, Jennifer L. Guler, Feifei Xiao

# Description
Copy number variants (CNVs) are prevalent in both diploid and haploid genomes. Studying CNVs in genomes is significantly advancing our knowledge of human disorders and disease susceptibility. Low-input data, including low-cell and single-cell sequencing, generally displays shallow and highly uneven read counts resulting from the whole genome amplification step that introduces amplification biases. Challenges also reside in selecting reference samples, which are used to provide baseline for defining copy number losses or gains. Traditionally, references are pre-specified from cells that are assumed to be normal or disease-free, which may bias CNV detection results if common CNVs are present among both tested and reference cells. Here we develop a new software and statistical framework for CNV detection in both haploid and diploid organisms with a novel strategy for defining references with low-input DNA sequencing data. The prominent advancement is the construction of a novel pseudo-reference to define biologically genuine references that results in effectively preserved common CNVs.


# Installation
GapCNV requires several R packages which requires manually installed: `modSaRa`,`FLCNA`,`GenomicRanges`, and `DNAcopy`. Then GapCNV can be installed by running:

```
install.packages("devtools")
library(devtools)
install_github("xuanxuanyu-bios/GapCNV")
```

# Running GapCNV
## Examples
Quality control and normalization for GC content and mappability
```
RC_QC <- Calling_rate_QC(RC=count.mat,ref_qc=ref_qc, 
                         cr.threshold=0.8,
                         GC.low=0.1,
                         GC.up=0.4,
                         mapp.threshold=0.9)
norm.data <- GC.MAP.normalization(RC=RC_QC$RC, ref_qc=RC_QC$ref_qc)
log2R     <- norm.data$log2Rdata
```
Construct pseudo-reference for each cell and perform data normalization
```
pseudo_ref<-GapCNV(count.mat=norm.data$RC_norm,
                 log2R=log2R,
                 ref=RC_QC$ref_qc,
                 nclust  = 2,
                 lambda  = 10,
                 cutoff  = 0.35)
```
Conduct CNV profiling using circular binary segmentation (CBS) method
```
CBS.res<-CBS.function(data.mat=pseudo_ref$log2R.norm.mat,
                  sample.name=NULL,alpha=0.01,
                  chr=as.vector(RC_QC$ref_qc@seqnames),
                  maploc=as.vector(RC_QC$ref_qc@ranges@start))
```
