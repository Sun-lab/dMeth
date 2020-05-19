
library(data.table)

# ------------------------------------------------------------
# read in probe information
# ------------------------------------------------------------

setwd("~/research/Deconvolution/data")

dat = fread("methylation_pure_ct_data.txt.gz")
datNew = fread("methylation_pure_ct_rmPC2_data.txt.gz")

dim(dat)
dim(datNew)

dat[1:2,1:3]
datNew[1:2,1:3]

dat = signif(dat, digits=4)
datNew = signif(datNew, digits=4)

dim(dat)
dim(datNew)

dat[1:2,1:3]
datNew[1:2,1:3]

fwrite(dat, file = "methylation_pure_ct_data_signif4.txt")
fwrite(datNew, file = "methylation_pure_ct_rmPC2_data_signif4.txt")

quit(save="no")


