
library(nnls)
library(glmnet)
library(data.table)

source("~/research/Deconvolution/R_batch1/_lib.R")

setwd("~/research/Deconvolution/data")

# ------------------------------------------------------------
# read in pure cell type data
# ------------------------------------------------------------

info = read.table("methylation_pure_ct_info.txt", as.is=TRUE, sep="\t",
                  header=TRUE, quote="")
dim(info)
info[1:2,]

dat = fread("methylation_pure_ct_rmPC2_data.txt")
dim(dat)
dat[1:2,1:5]

sam = read.table("methylation_pure_ct_sample.txt", as.is=TRUE,
                 sep="\t", header=TRUE)
dim(sam)
sam[1:2,]

table(names(dat) == sam$id)

dat = data.matrix(dat)
table(sam$label)
rownames(dat) = info$ID

# ------------------------------------------------------------
# read DNA methylation data
# ------------------------------------------------------------

setwd("~/research/TCGA/COAD/")

datM = fread(file = "_data2/methylation_betaValue.txt", header=TRUE)
dim(datM)
datM[1:2, 1:5]

infoM = read.table(file = "_data2/methylation_info.txt", sep = "\t",
                   header = TRUE, as.is = TRUE, quote="")
dim(infoM)
infoM[1:2, 1:5]
table(infoM$CHR)

table(datM$id == infoM$Name)

datM = data.matrix(datM[,-1])
rownames(datM) = infoM$Name
dim(datM)
datM[1:2, 1:5]

# ------------------------------------------------------------
# take intersection
# ------------------------------------------------------------

table(info$ID == info$Name)
cpgs = intersect(info$ID, infoM$Name)
length(cpgs)

mat1  = match(cpgs, infoM$Name)
datM  = datM[mat1,]
infoM = infoM[mat1,]

dim(datM)
datM[1:2, 1:5]

dim(infoM)
infoM[1:2, 1:5]

mat2 = match(cpgs, info$ID)
dat  = dat[mat2,]
info = info[mat2,]

dim(dat)
dat[1:2, 1:5]

dim(info)
info[1:2, 1:5]

# ------------------------------------------------------------
# read DNA methylation sample information
# ------------------------------------------------------------

setwd("~/research/TCGA/COAD/_data2/")

coad = read.table("patient_coad_short_table_nMut.txt", sep="\t",
as.is=TRUE, header=TRUE)
dim(coad)
coad[1:2,]

ff0    = "patient_coad_M_info_hyperMeth.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is = TRUE)
dim(emInfo)
emInfo[1, ]

table(emInfo$bcr_patient_barcode %in% coad$bcr_patient_barcode)
mat1 = match(emInfo$bcr_patient_barcode, coad$bcr_patient_barcode)
emInfo$nMut = coad$nMut[mat1]

table(colnames(datM) == emInfo$patient_id)

# ------------------------------------------------------------
# read in probes to be used
# ------------------------------------------------------------

setwd("~/research/Deconvolution/data")

load("methylation_p2use.RData")
lapply(p2use, dim)

nms = names(p2use)
for(nm1 in nms){
  p1 = p2use[[nm1]]
  p2use[[nm1]] = p1[which(p1$id %in% cpgs),]
}

lapply(p2use, dim)

genes = NULL
for(lb1 in names(p2use)){
  genes = c(genes, p2use[[lb1]]$id)
}

length(genes)
length(unique(genes))
genes = unique(genes)

table(genes %in% rownames(datM))
table(genes %in% rownames(dat))

X  = dat[match(genes, rownames(dat)),]
dim(X)
colnames(X)

dim(sam)
table(sam$id == colnames(X))
table(sam$label)

# ------------------------------------------------------------
# extract methylation data from tumor samples
# ------------------------------------------------------------

ys = datM[match(rownames(X), rownames(datM)),]
dim(ys)
ys[1:2,1:5]

# ------------------------------------------------------------
# non-negative least square
# ------------------------------------------------------------

cellTypes = unique(sam$label)

rhos1 = rhos2 = matrix(NA, nrow=ncol(ys), ncol=length(cellTypes))
nIt1  = nIt2  = nobs = rep(NA, ncol(ys))

for(j in 1:ncol(ys)){
  if(j %% 50 == 0){ cat(j, date(), "\n") }
  y    = ys[,j]
  
  wj   = wls(y, X, cellTypes, sam)
  rhos1[j,] = wj$b2
  nIt1[j]   = wj$nIt
  
  wj   = wlsEM(y, X, cellTypes, sam)
  rhos2[j,] = wj$b2
  nIt2[j]   = wj$nIt
  nobs[j]   = wj$nobs
}

summary(nIt1)
summary(nIt2)
summary(nobs)

rownames(rhos1) = rownames(rhos2) = colnames(ys)
colnames(rhos1) = colnames(rhos2) = cellTypes

cor(rhos1, emInfo$nMut, use="pair")
apply(rhos1, 2, cor.test, y=emInfo$nMut)

cor(rhos1, emInfo$abs_purity, use="pair")
cor(rhos1, emInfo$hyperMeth, use="pair")

setwd("~/research/Deconvolution/figures")

pdf("wls_vs_wlsEM.pdf", width=4, height=4)
par(mar=c(5,4,1,1), bty="n")
plot(rhos1, rhos2, xlab="wls", ylab="wlsEM")
dev.off()

pdf("heatmap_wls.pdf", width=4, height=6)
heatmap(rhos1, margins=c(3,1), labRow="")
dev.off()

pdf("heatmap_wlsEM.pdf", width=4, height=6)
heatmap(rhos2, margins=c(3,1), labRow="")
dev.off()

# ------------------------------------------------------------
# non-negative least square, with purity information
# ------------------------------------------------------------

rhos1P = rhos2P = matrix(NA, nrow=ncol(ys), ncol=length(cellTypes))
nIt1   = nIt2   = nobs = rep(NA, ncol(ys))

summary(emInfo$abs_purity)
table(emInfo$abs_purity > 0.9)

for(j in 1:ncol(ys)){
  if(j %% 50 == 0){ cat(j, date(), "\n") }
  y      = ys[,j]
  purity = emInfo$abs_purity[j]
  
  if(is.na(purity)){ next }
  if(purity > 0.9){ next }

  wj   = wls(y, X, cellTypes, sam, total=1-purity)
  rhos1P[j,] = wj$b2
  nIt1[j]    = wj$nIt
  
  wj   = wlsEM(y, X, cellTypes, sam, total=1-purity)
  rhos2P[j,] = wj$b2
  nIt2[j]    = wj$nIt
  nobs[j]    = wj$nobs
}

summary(nIt1)
summary(nIt2)
summary(nobs)

rownames(rhos1P) = rownames(rhos2P) = colnames(ys)
colnames(rhos1P) = colnames(rhos2P) = cellTypes

cor(rhos1P, emInfo$nMut, use="pair")
apply(rhos1P, 2, cor.test, y=emInfo$nMut)

cor(rhos1P, emInfo$abs_purity, use="pair")
cor(rhos1P, emInfo$hyperMeth, use="pair")

w2kp = which(emInfo$abs_purity <= 0.9)

setwd("~/research/Deconvolution/figures")

pdf("wls_vs_wlsEM_given_purity.pdf", width=4, height=4)
par(mar=c(5,4,1,1), bty="n")
plot(rhos1P[w2kp], rhos2P[w2kp], xlab="wls", ylab="wlsEM")
dev.off()

pdf("heatmap_wls_given_purity.pdf", width=4, height=6)
heatmap(rhos1P[w2kp,], margins=c(3,1), labRow="")
dev.off()

pdf("heatmap_wlsEM_given_purity.pdf", width=4, height=6)
heatmap(rhos2P[w2kp,], margins=c(3,1), labRow="")
dev.off()

# ------------------------------------------------------------
# elastic net
# ------------------------------------------------------------

rhos3 = rhos3P = matrix(NA, nrow=ncol(ys), ncol=length(cellTypes))
nIt3  = nIt3P  = rep(NA, ncol(ys))

for(j in 1:ncol(ys)){
  if(j %% 10 == 0){ cat(j, date(), "\n") }
  y  = ys[,j]
  wj = enet(y, X, cellTypes, sam, alpha=0.2)

  rhos3[j,] = wj$b2
  nIt3[j]   = wj$nIt
  
  purity = emInfo$abs_purity[j]
  
  if(is.na(purity)){ next }
  if(purity > 0.9){ next }
  
  wj = enet(y, X, cellTypes, sam, alpha=0.2, total=1-purity)
  rhos3P[j,] = wj$b2
  nIt3P[j]   = wj$nIt
}

summary(nIt3)
summary(nIt3P)

rownames(rhos3) = rownames(rhos3P) = colnames(ys)
colnames(rhos3) = colnames(rhos3P) = cellTypes

cor(rhos3, emInfo$nMut, use="pair")
cor(rhos3, emInfo$abs_purity, use="pair")
cor(rhos3, emInfo$hyperMeth, use="pair")

cor(rhos3P, emInfo$nMut, use="pair")
cor(rhos3P, emInfo$abs_purity, use="pair")
cor(rhos3P, emInfo$hyperMeth, use="pair")

pdf("heatmap_enet.pdf", width=4, height=6)
heatmap(rhos3, margins=c(3,1), labRow="")
dev.off()

pdf("heatmap_enet_given_purity.pdf", width=4, height=6)
heatmap(rhos3P[w2kp,], margins=c(3,1), labRow="")
dev.off()

# ------------------------------------------------------------
# save results
# ------------------------------------------------------------

save(rhos1, rhos2, rhos3, file="../data/methylation_deconv.Rdata")
save(rhos1P, rhos2P, rhos3P, file="../data/methylation_deconv_given_purity.Rdata")

quit(save="no")


