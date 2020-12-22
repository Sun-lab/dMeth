# Pipeline for the BLUEPRINT data

This is the pipelines to use the BLUEPRINT data (Chen et al. (2016)) to construct *in silico* mixture using cell type-specific DNA methylation and gene expression data in three human immune cell types (Monocytes, Neutrophils and naive CD4 T cells) from 197 individuals. Details of experimental design and results are provided in the supplementary materials of our paper (Zhang et al. (2021)). This pipeline includes:

* step0-GetData.R: preprocess the methylation data matrix.

* step1-estMeanVar.R: estimate the mean and variance from methylation and remove probes with extremely large variance. Then separate the data and obtain ```dat.gen``` to generate mixture and ```dat.est``` to construct the reference data.

* step2-filterProbes.R:  select probes based on ```dat.est``` and save in 'probes.RData'

* step3-RunMethy.R: generate methylation mixture and run cell type deconvolution. 

* step4-runExpr.R: generate gene expression mixture data and then perform CIBERSORT and compare the result with results from methylation data.

# Reference

Zhang et al. 2021, EMeth: An EM algorithm for cell type decomposition based on DNA methylation data.

Chen et al. (2016) Genetic drivers of epigenetic and transcriptional variation in human immune cells. Cell, 167(5), 1398-1414.

