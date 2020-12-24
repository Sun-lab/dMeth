# TCGA Real Data Analysis

This pipelines are codes to run real data analysis using TCGA data. The files include:

- The procedure for preparing gene expression data are included in 1\_eData\_\{CancerType\}.Rmd files.  The processed gene expression data are then uploaded into CIBERSORTx website to obtain estimates of cell type proportions for each bulk tumor sample. 

- 2\_run\_COAD.R and 2\_run\_another\_cancer_type.R:  estimates of cell type proportions for each cancer type. The code for COAD is a bit different from other cancer types because it has its own input and output file formats. The other file 2\_run_another_cancer_type.R can be used to run multiple cancer types using the shell script  2_run_another_cancer_type.sh. These files include the following steps.

  - read in pure cell type data: data are stored in ../cell_type_specific_reference/data, include

    - methylation_pure_ct_info.txt.gz
    - methylation_pure_ct_rmPC2_data_signif4.txt.gz
    - methylation_pure_ct_sample.txt

  - read in DNA methylation data: the source files are

    - methylation_betaValue.txt
    - methylation_info.txt

  - take intersection of CpG probes between purified data and bulk tumor samples.

  - read in probes to be used: ref_966probes.RData, located in

    ../cell_type_specific_reference/data, this is also used in the Simulation Study.

  - Compute and compare estimation results from DNA methylation data

- \{CancerType\}_ct_prop_comparison.R: compares estimates of cell type proportions for each cancer type from Methylation and gene expression data. And generate corresponding plots.

