# Pipeline for the BLUEPRINT data

This is the pipelines to use the BLUEPRINT data (Chen et al. (2016)) to construct *in silico* mixture using cell type-specific DNA methylation and gene expression data in three human immune cell types (Monocytes, Neutrophils and naive CD4 T cells) from 197 individuals. Details of experimental design and results are provided in the supplementary materials of our paper (Zhang et al. (2021)). This pipeline includes:

* step0-GetData.R: preprocess the methylation data matrix.
  * Source data can be downloaded from [ftp://ftp.ebi.ac.uk/pub/databases/blueprint/blueprint_Epivar/Pheno_Matrix/](ftp://ftp.ebi.ac.uk/pub/databases/blueprint/blueprint_Epivar/Pheno_Matrix/)
  * Add normal noise to beta value of methylation data
* step1-estMeanVar.R: estimate the mean and variance from methylation and remove probes with extremely large variance. Then separate the data and obtain ```dat.gen``` to generate mixture and ```dat.est``` to construct the reference data.
  * load perturbed data (with normal noise)
  * filter probe to keep those with small variance
  * estimate the mean and variance from methylation
  * seperate data and obtain `dat.gen` to generate mixture and `dat.est` to construct the reference data.
* step2-filterProbes.R:  select probes based on ```dat.est``` and save in 'probes.RData' based on differential methylation. P-value was set to be 1e-6.
* step3-RunMethy.R: 
  * generate methylation mixture 
  * run cell type deconvolution
  * results are saved in
    * cormat_\{perturblevel\}.RData, the correlation between true proportion and estimated proportion
    * err_\{perturblevel\}.RData: the MSE for each cell type
    * samp_est.RData: the sample used to estimate cell type proportion
    * samp_gen.RData: the sample used to generate mixture and filter probes
    * truerho.RData: the true proportion used when generating mixtures
    * methyrho.RData: the estimation result of cell type proportion
* step4-runExpr.R: 
  * Load gene expression data, add same level of noise compared to the methylation data.
  * filter gene with differential gene expression level and p value equals 1e-6.
  * seperate cell type specific gene expression data into `samp_gen` and `samp_est` and choose common samples in both expression data and methylation data.
  * generate pseudo mixtures of gene expression data and run CIBERSORT source code to estimate the cell type proportion. CIBERSORT source code can be downloaded from the website [https://cibersort.stanford.edu](https://cibersort.stanford.edu). 
  * compare the result of two estimations.

# Reference

Zhang et al. 2021, EMeth: An EM algorithm for cell type decomposition based on DNA methylation data.

Chen et al. (2016) Genetic drivers of epigenetic and transcriptional variation in human immune cells. Cell, 167(5), 1398-1414.

