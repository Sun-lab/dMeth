# Deconvolute cell type composition using DNA methylation data. 

This repository includes the pipelines for simulation and real data analysis for our paper (Zhang et al. (2021)).  It includes four subdirectories, each representing a different pipeline. 

To test the codes or reproduce part of our analysis, please first clone this repository to your local machine. Then the source file in these folders can be directly executed in R. These pipelines use our EMeth package, which requires R version 3.6.0 or higher and can be installed by the following codes in R.

 ```
    library("devtools");
    install_github("Hanyuz1996/EMeth")
 ```
The subdirectories includes codes for different steps. 

1. Blueprint_pipeline: source code for the *in silico* mixture study using the BLUEPRINT data. 

2. Simulation_pipeline: source code for the simulation study, where we simulated mixture based on cell type-specific reference data. 

3. TCGA_pipeline: source code for the TCGA study, comparing cell type deconvolution from gene expression data and DNA methylation data.

4. cell_type_specific_reference: source data and codes for preparation for cell type specific reference data. This reference data is used by both  simulation study and TCGA study. 


Details of each pipeline can be found in README.md in each subdirectory.

# Reference

Zhang et al. 2021, EMeth: An EM algorithm for cell type decomposition based on DNA methylation data.

