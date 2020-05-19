
# the raw data of DNA metylation is too large to be kept in gitHub
# here is the local path for DNA methylation data
path.data = "~/research/TCGA/COAD/_data2"

library(data.table)
library(stringr)

# ------------------------------------------------------------
# read in pure cell type data
# ------------------------------------------------------------

path.ref = "../cell_type_specific_reference/data"

info = fread(file.path(path.ref, "methylation_pure_ct_info.txt.gz"))
dim(info)
info[1:2,]

dat = fread(file.path(path.ref, 
                      "methylation_pure_ct_rmPC2_data_signif4.txt.gz"))
dim(dat)
dat[1:2,1:5]

sam = fread(file.path(path.ref, "methylation_pure_ct_sample.txt"))
dim(sam)
sam[1:2,]

table(names(dat) == sam$id)

dat = data.matrix(dat)
table(sam$label)
rownames(dat) = info$ID

# ------------------------------------------------------------
# read DNA methylation data
# ------------------------------------------------------------

datM = fread(file.path(path.data, "methylation_betaValue.txt"))
dim(datM)
datM[1:2, 1:5]

infoM = fread(file.path(path.ref, "methylation_info.txt"))
dim(infoM)
names(infoM)
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

ff0    = "../TCGA_results/clinical_data/patient_coad_M_info_hyperMeth.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is = TRUE)
dim(emInfo)
emInfo[1, ]

# ------------------------------------------------------------
# read in probes to be used
# ------------------------------------------------------------

load(file.path(path.ref, "ref_966probes.RData"))
ls()
genes = probe2use
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
cellTypes = unique(sam$label)

# ------------------------------------------------------------
# extract methylation data from tumor samples
# ------------------------------------------------------------

ys = datM[match(rownames(X), rownames(datM)),]
dim(ys)
ys[1:2,1:5]

ys_na    = which(apply(is.na(ys),2,any))
eta_abs  = emInfo$abs_purity
eta_abs_na = which(is.na(eta_abs))

any.na = union(ys_na,eta_abs_na)
any.na

ys = ys[,-any.na]
eta_abs = eta_abs[-any.na]

dim(ys)

#-------------------------------------------------------------
# Estimate Mean Matrix mu
#-------------------------------------------------------------

mu = matrix(NA, nrow = dim(X)[1], ncol = length(cellTypes))
s2 = matrix(NA, nrow = dim(X)[1], ncol = length(cellTypes))

row.names(mu) = row.names(s2) = rownames(X)
colnames(mu)  = colnames(s2)  = cellTypes

for(ct in cellTypes){
  sam.ct = unlist(sam[which(sam[,2]==ct),1])
  dat.ct = X[,sam.ct]
  mu[,ct] = rowMeans(dat.ct,na.rm=TRUE)
  s2[,ct] = apply(dat.ct,1,sd,na.rm = TRUE)^2
}

#-------------------------------------------------------------
# Read Estimation from Expression Data
#-------------------------------------------------------------

fnm = '../TCGA_results/deconv_results/COAD_composition_cibersortx.txt'
est_expr = fread(fnm)
dim(est_expr)
est_expr[1:2,]

samname  = str_replace(est_expr$Mixture, "^X", "")
length(samname)
samname[1:5]

com_sample = intersect(samname,colnames(ys))
length(com_sample)

w2kp = which(samname %in% com_sample)
est_expr = data.matrix(est_expr[w2kp,-1])
rownames(est_expr) = samname[w2kp]
dim(est_expr)
est_expr[1:2,]

deconv_expr = matrix(NA, nrow = nrow(est_expr),ncol = length(cellTypes))
colnames(deconv_expr) = cellTypes 
rownames(deconv_expr) = rownames(est_expr)
colnames(est_expr)

other = rowSums(est_expr[,c(10,19,20,21)])
deconv_expr[,"B"] = rowSums(est_expr[,1:3])/0.4
deconv_expr[,"CD4T"] = rowSums(est_expr[,5:8])/0.4
deconv_expr[,"CD8T"] = as.matrix(est_expr[,4])/0.4
deconv_expr[,"Treg"] = as.matrix(est_expr[,9])/0.4
deconv_expr[,"NK"] = rowSums(est_expr[,11:12])/0.42
deconv_expr[,"Monocyte"] = rowSums(est_expr[,13:18])/1.40
deconv_expr[,"Neutrophil"] = as.matrix(est_expr[,22])/0.15
deconv_expr = deconv_expr / rowSums(deconv_expr)

print(dim(deconv_expr))
deconv_expr[1:2,]






