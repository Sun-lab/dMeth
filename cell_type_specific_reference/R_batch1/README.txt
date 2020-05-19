Steps to process the methylation data. Appologize that the path was set for the file system in my laptop and was not updated according to the GitHub. Most of the data in folder "~/research/Deconvolution/data" were copied to folder ../data, except some files larger than 100Mb after gzip, since GitHub does not allow a single file larger than 100Mb. Most of data in folder "~/research/Deconvolution/original_data" are too large, so only part of the data were copied to folder ../original_data so that they can be used to test the code. 

-----------------------------------------------------------------------------
step1_check_ENA_samples.R
-----------------------------------------------------------------------------

	Query the ENA database for cell type-specific DNA methylation data and saved the query results in folder _ENA_by_cell_type. Then step1_check_ENA_samples.R checked these query results to select apporpriate ones for the next step

-----------------------------------------------------------------------------
step2_get_data.R
-----------------------------------------------------------------------------

	To select a set of CpG pobes that are in genome build 37, with CHR information and does not have a SNP within 10bp. Then extract data of these CpGs across all the studies. 
	
    > table(i450$Genome_Build)

               36     37 
        65   3091 482421 
    > table(i450$CHR)

        Y     X     1     2     3     4     5     6     7     8     9    10    11 
      416 11232 46857 34810 25159 20464 24327 36611 30017 20950  9861 24388 28794 
       12    13    14    15    16    17    18    19    20    21    22       
    24539 12285 15078 15259 21969 27879  5922 25521 10379  4243  8552    65 
    > 
    > table(i450$Methyl27_Loci)

             TRUE 
    459500  25978 
    > table(i450$Probe_SNPs=="")

     FALSE   TRUE 
     59892 425685 
    > table(i450$Probe_SNPs_10=="")

     FALSE   TRUE 
     36535 449042 
    > 
    > w2use = (i450$Genome_Build=="37" & i450$CHR != "")
    > w2use = w2use & (i450$Probe_SNPs=="" & i450$Probe_SNPs_10=="")
    > 
    > table(w2use)
    w2use
     FALSE   TRUE 
     92176 393401 

-----------------------------------------------------------------------------
step3_process_methylation_data.R
-----------------------------------------------------------------------------

1. After initial checking by PCA (plots methylation_PCA_eigen_values.pdf, methylation_PCA_PC1_vs_PC2.pdf etc.), remove one study Vento of 9 samples, further remove two cell types: "Eosinophils", "Granulocytes". 


2. remove batch effect at m-value scale, using the data from Reinius et al. as reference. The reason we used Reinius et al. as reference is that it includes the data from different cell types. 
2.1 quantile normalization on the data from Reinius et al.   
2.2 align the data from other study to the quantile normalized data from Reinius et al. 

Generate the set of plots methylation_PCA_after_nm_*. Obtain the following set of files:

methylation_pure_ct_data.txt:		Pure cell type data file
methylation_pure_ct_sample.txt:	Pure cell type sample file
methylation_pure_ct_info.txt:		Pure cell type information file, of the selected probes

3. Joint PCA, and PC2 is related with studies. Regress out PC2

obtain the following files: 

methylation_pure_ct_rmPC2_data.txt.txt:		Pure cell type data file

-----------------------------------------------------------------------------
step3.0_reduce_file_size.R
-----------------------------------------------------------------------------

reduce the file size by keeping 4 digits


