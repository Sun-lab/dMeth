# TCGA Real Data Analysis

This pipelines are codes to run real data analysis using TCGA data. The files include:

- runCOAD.R, runLUAD.R, runLUSC.R, runSKCM.R: preprocess methylation data for each cancer type and run EMeth and benchmark algorithms and compare with previously estimated result from CIBERSORTX. 

- Cell type fraction estimated by cibersortx using gene expression data were saved in folder _cibersortx_results.
