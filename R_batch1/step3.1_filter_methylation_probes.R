
source("~/research/Deconvolution/R_batch1/_lib.R")

setwd("~/research/Deconvolution/data")

library(multcomp)
library(data.table)

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------

info = read.table("methylation_pure_ct_info.txt", as.is=TRUE, sep="\t",
                  header=TRUE, quote="")
dim(info)
info[1:2,]

dat = fread("methylation_pure_ct_rmPC2_data.txt")
dim(dat)
dat[1:2,1:5]

dat0 = fread("methylation_pure_ct_data.txt")
dim(dat0)
dat0[1:2,1:5]

sam = read.table("methylation_pure_ct_sample.txt", as.is=TRUE,
    sep="\t", header=TRUE)
dim(sam)
sam[1:2,]

table(names(dat) == sam$id)

dat  = data.matrix(dat)
dat0 = data.matrix(dat0)
table(sam$label)

rownames(dat)  = info$ID
rownames(dat0) = info$ID

cr1 = cor(dat, dat0, use="pair")
summary(diag(cr1))
summary(cr1[upper.tri(cr1)])

pdf("../figures/before_vs_after_rmPC2.pdf", width=4, height=4)
par(mar=c(5,4,1,1), bty="n")
smoothScatter(as.numeric(dat0), as.numeric(dat), xlab="before", ylab="after")
abline(0,1)
dev.off()

# ------------------------------------------------------------
# define the cell types to separate
# ------------------------------------------------------------

panT = c("CD4T", "CD8T", "Treg")
p2use = levels = list()

level1 = list()
level1[["lymphoid"]] = c("B", "CD4T", "CD8T", "Treg", "NK")
level1[["myeloid"]]  = c("Monocyte", "Neutrophil")

levels[[1]] = level1

level1 = list()
level1[["B_NK"]] = c("B", "NK")
level1[["T"]]    = panT

levels[[2]] = level1

level1 = list()
level1[["B"]]  = "B"
level1[["NK"]] = "NK"

levels[[3]] = level1

level1 = list()
level1[["CD4T_Tregs"]] = c("CD4T", "Tregs")
level1[["CD8T"]] = "CD8T"

levels[[4]] = level1

level1 = list()
level1[["CD4T"]] = "CD4T"
level1[["Tregs"]]   = "Tregs"

levels[[5]] = level1

level1 = list()

level1[["Neutrophil"]] = "Neutrophil"
level1[["Monocyte"]]   = "Monocyte"

levels[[6]] = level1

uct = unique(sam$label)
uct

for(k in 1:length(uct)){
  ct1 = uct[k]
  ct2 = paste("Non", ct1)
  
  level1 = list()

  level1[[ct1]] = ct1
  level1[[ct2]] = setdiff(uct, ct1)
  
  levels[[6+k]] = level1
}

length(levels)

path  = "../figures_methylation"
dataType = "methylation"

for(kk in 1:length(levels)){
  cat(kk, date(), "\n")
  level1 = levels[[kk]]
  nms    = names(level1)
  nm1    = sprintf("%s_vs_%s", nms[1], nms[2])
  
  p2use[[nm1]] = probeSelect(level1, sam, dat, path=path, dataType=dataType)
}

lapply(p2use, dim)

setwd("~/research/Deconvolution/data")

save(p2use, file="methylation_p2use.RData")

quit(save="no")


