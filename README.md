# Deconvolute cell type composition using DNA methylation data. 

This repository includes the source code and part of the original dataset we used for the paper _To insert the paper link once the paper is online_. It includes four folders subdirectories, each representing a different pipeline for the paper. 

To test the codes for repeating the study, first clone this repository to your local machine. And the source file in these folders can be directly executed in R. The code require installation of EMeth package, which requires R version 3.6.0 or higher by the following codes in R.

 ```
    library("devtools");
    install_github("Hanyuz1996/EMeth")
 ```
The subdirectories includes codes for different steps. 

1. cell_type_specific_reference: source data and codes for preparation for cell type specific reference data. This reference data is used both for CombinedStudy and TCGA study.
2. CombinedStudy_pipeline: source code for the simulation study.
3. TCGA_pipeline: source code for the TCGA study, comparing cell type deconvolution from gene expression data and DNA methylation data.
4. Blueprint_pipeline: source code for the blueprint data study.

Details of each pipeline can be found in README.md in each subdirectories.