
library(data.table)
library(ggplot2)
# library(cowplot)

fun1 <- function(v, f){
  tapply(v, f, mean, na.rm=TRUE)
}

# ------------------------------------------------------------
# read in probe information
# ------------------------------------------------------------

setwd("~/research/Deconvolution/data")

info = fread("i450_filtered_probes.txt")
dim(info)
info[1:2,]

info = info[,-c(3:6,10,14)]
dim(info)
info[1:2,]

g1 = grep("FOXP3", info$UCSC_RefGene_Name)
info[g1,]

Foxp3 = info$ID[g1]
Foxp3

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------

sam.Reinius = fread("Reinius_methy_samples.txt")
dat.Reinius = fread("Reinius_methy_data.txt")

dim(dat.Reinius)
dat.Reinius[1:3,1:5]

dim(sam.Reinius)
sam.Reinius[1:3,]
table(sam.Reinius$label)

nNA.R = rowSums(is.na(dat.Reinius))
table(nNA.R)

dat.Reinius = dat.Reinius[which(nNA.R <= 2),]
dim(dat.Reinius)
dat.Reinius[1:3,1:5]

table(sam.Reinius$label)
ctMeans = t(apply(dat.Reinius[,-1,with=FALSE], 1, fun1, f=sam.Reinius$label))
dim(ctMeans)
ctMeans[1:2,]

ww1 = which(is.na(dat.Reinius), arr.ind=TRUE)
dim(ww1)
ww1[1:5,]

for(i in 1:nrow(ww1)){
  if(i %% 5000 == 0){ cat(i, date(), "\n") }
  rid = ww1[i,1]
  cid = ww1[i,2]
  
  mati = which(colnames(ctMeans) == sam.Reinius$label[cid-1])
  dat.Reinius[[cid]][rid] = ctMeans[rid,mati]
}

table(Foxp3 %in% dat.Reinius$id)
mat1 = match(dat.Reinius$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# CD4T_CD8T_E-GEOD-71957
# ------------------------------------------------------------

dat.Limbach = fread("CD4T_CD8T_E-GEOD-71957_methy_data.txt")
dim(dat.Limbach)
dat.Limbach[1:2,1:5]

sam.Limbach = fread("CD4T_CD8T_E-GEOD-71957_methy_sample.txt")
dim(sam.Limbach)
sam.Limbach[1:2,]
table(sam.Limbach$cell_type)

nNA.R = rowSums(is.na(dat.Limbach))
table(nNA.R)

dim(dat.Limbach)
dat.Limbach = dat.Limbach[which(nNA.R <= 2),]
dim(dat.Limbach)
dat.Limbach[1:3,1:5]

nNA.R = rowSums(is.na(dat.Limbach))
table(nNA.R)

table(Foxp3 %in% dat.Limbach$id)
mat1 = match(dat.Limbach$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# mono_CD4T_E-GEOD-56047, read in sample information
# ------------------------------------------------------------

sam.Reynolds = fread("mono_CD4T_E-GEOD-56047_methy_sample.txt")
dim(sam.Reynolds)
sam.Reynolds[1:2,]

apply(sam.Reynolds[,5:9], 2, summary)

length(unique(sam.Reynolds$altID))
table(table(sam.Reynolds$altID))

table(sam.Reynolds$well_id)
table(sam.Reynolds$Sample_source_name)
table(sam.Reynolds$Sample_source_name, sam.Reynolds$well_id)

ww1 = which(sam.Reynolds$Sample_source_name == "CD4+ cell")
sam.Reynolds.CD4  = sam.Reynolds[ww1]
sam.Reynolds.mono = sam.Reynolds[-ww1]

apply(sam.Reynolds.CD4[,5:9], 2, summary)
apply(sam.Reynolds.mono[,5:9], 2, summary)

max.CD4 = apply(sam.Reynolds.CD4[,5:9], 1, max, na.rm = TRUE)
summary(max.CD4)
table(max.CD4 < 0.2)

max.CD14 = apply(sam.Reynolds.mono[,5:9], 1, max, na.rm = TRUE)
summary(max.CD14)
table(max.CD14 < 0.15)

sam.Reynolds.CD4 = sam.Reynolds.CD4[which(max.CD4 < 0.2),]
dim(sam.Reynolds.CD4)
sam.Reynolds.CD4[1:2,]

sam.Reynolds.mono = sam.Reynolds.mono[which(max.CD14 < 0.15),]
dim(sam.Reynolds.mono)
sam.Reynolds.mono[1:2,]

# ------------------------------------------------------------
# CD4T data
# ------------------------------------------------------------

dat.Reynolds.CD4 = fread("CD4T_E-GEOD-56047_methy_data.txt")
dim(dat.Reynolds.CD4)
dat.Reynolds.CD4[1:2,1:5]

mat1 = match(sam.Reynolds.CD4$Assay_Name, names(dat.Reynolds.CD4))
dat.Reynolds.CD4 = dat.Reynolds.CD4[,c(1,mat1),with=FALSE]
dim(dat.Reynolds.CD4)
dat.Reynolds.CD4[1:2,1:5]

table(names(dat.Reynolds.CD4)[-1] == sam.Reynolds.CD4$Assay_Name)

summary(dat.Reynolds.CD4[[2]])

for(k in 2:ncol(dat.Reynolds.CD4)){
  nmk = names(dat.Reynolds.CD4)[k]
  dat.Reynolds.CD4[[nmk]] = exp(dat.Reynolds.CD4[[nmk]])/(1 + exp(dat.Reynolds.CD4[[nmk]]))
}
summary(dat.Reynolds.CD4[[2]])

nNA.R = rowSums(is.na(dat.Reynolds.CD4))
table(nNA.R)

names(dat.Reynolds.CD4)[1] = "id"
table(Foxp3 %in% dat.Reynolds.CD4$id)
mat1 = match(dat.Reynolds.CD4$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# mono data
# ------------------------------------------------------------

dat.Reynolds.mono = fread("CD14T_E-GEOD-56047_methy_data.txt")
dim(dat.Reynolds.mono)
dat.Reynolds.mono[1:2,1:5]

mat1 = match(sam.Reynolds.mono$Assay_Name, names(dat.Reynolds.mono))
dat.Reynolds.mono = dat.Reynolds.mono[,c(1,mat1),with=FALSE]
dim(dat.Reynolds.mono)
dat.Reynolds.mono[1:2,1:5]

table(names(dat.Reynolds.mono)[-1] == sam.Reynolds.mono$Assay_Name)

summary(dat.Reynolds.mono[[2]])
for(k in 2:ncol(dat.Reynolds.mono)){
  nmk = names(dat.Reynolds.mono)[k]
  dat.Reynolds.mono[[nmk]] = exp(dat.Reynolds.mono[[nmk]])/(1 + exp(dat.Reynolds.mono[[nmk]]))
}
summary(dat.Reynolds.mono[[2]])

nNA.R = rowSums(is.na(dat.Reynolds.mono))
table(nNA.R)

names(dat.Reynolds.mono)[1] = "id"
table(Foxp3 %in% dat.Reynolds.mono$id)
mat1 = match(dat.Reynolds.mono$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# mono_DC_MAC_E-GEOD-71837
# ------------------------------------------------------------

dat.Vento = fread("mono_DC_MAC_E-GEOD-71837_methy_data.txt")
dim(dat.Vento)
dat.Vento[1:2,1:5]

sam.Vento = fread("mono_DC_MAC_E-GEOD-71837_methy_sample.txt")
dim(sam.Vento)
sam.Vento[1:2,]
table(sam.Vento$cell_type)

nNA.R = rowSums(is.na(dat.Vento))
table(nNA.R)

dim(dat.Vento)
dat.Vento = dat.Vento[which(nNA.R <= 1),]
dim(dat.Vento)
dat.Vento[1:3,1:5]

table(sam.Vento$cell_type)
ctMeans = t(apply(dat.Vento[,-1,with=FALSE], 1, fun1, f=sam.Vento$cell_type))
dim(ctMeans)
ctMeans[1:2,]

ww1 = which(is.na(dat.Vento), arr.ind=TRUE)
dim(ww1)
ww1[1:5,]

for(i in 1:nrow(ww1)){
  rid = ww1[i,1]
  cid = ww1[i,2]
  
  mati = which(colnames(ctMeans) == sam.Vento$cell_type[cid-1])
  dat.Vento[[cid]][rid] = ctMeans[rid,mati]
}

nNA.R = rowSums(is.na(dat.Vento))
table(nNA.R)

table(Foxp3 %in% dat.Vento$id)
mat1 = match(dat.Vento$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# neut_E-GEOD-65097_methy_data.txt
# ------------------------------------------------------------

dat.Coit = fread("neut_E-GEOD-65097_methy_data.txt")
dim(dat.Coit)
dat.Coit[1:2,1:5]

sam.Coit = fread("neut_E-GEOD-65097_methy_sample.txt")
dim(sam.Coit)
sam.Coit[1:2,]
table(sam.Coit$cell_type)

nNA.R = rowSums(is.na(dat.Coit))
table(nNA.R)

dim(dat.Coit)
dat.Coit = dat.Coit[which(nNA.R <= 3),]
dim(dat.Coit)
dat.Coit[1:3,1:5]

table(sam.Coit$cell_type)
ctMeans = rowMeans(dat.Coit[,-1,with=FALSE], na.rm=TRUE)
length(ctMeans)
ctMeans[1:5]

ww1 = which(is.na(dat.Coit), arr.ind=TRUE)
dim(ww1)
ww1[1:5,]

for(i in 1:nrow(ww1)){
  rid = ww1[i,1]
  cid = ww1[i,2]
  
  dat.Coit[[cid]][rid] = ctMeans[rid]
}

nNA.R = rowSums(is.na(dat.Coit))
table(nNA.R)

table(Foxp3 %in% dat.Coit$id)
mat1 = match(dat.Coit$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# NK_E-GEOD-66562_methy_data.txt
# there is no informtion on chrX or chrY
# ------------------------------------------------------------

dat.Schlums = fread("NK_E-GEOD-66562_methy_data.txt")
dim(dat.Schlums)
dat.Schlums[1:2,1:5]

sam.Schlums = fread("NK_E-GEOD-66562_methy_sample.txt")
dim(sam.Schlums)
sam.Schlums[1:2,]
table(sam.Schlums$cell_type)

nNA.R = rowSums(is.na(dat.Schlums))
table(nNA.R)

table(Foxp3 %in% dat.Schlums$id)
mat1 = match(dat.Schlums$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# Treg_E-GEOD-49667_methy_data.txt
# ------------------------------------------------------------

dat.Zhang = fread("Treg_E-GEOD-49667_methy_data.txt")
dim(dat.Zhang)
dat.Zhang[1:2,1:5]

sam.Zhang = fread("Treg_E-GEOD-49667_methy_sample.txt")
dim(sam.Zhang)
sam.Zhang[1:2,]

table(sam.Zhang$cell_type)

nNA.R = rowSums(is.na(dat.Zhang))
table(nNA.R)

table(Foxp3 %in% dat.Zhang$id)
mat1 = match(dat.Zhang$id, info$ID)
table(info$CHR[mat1])

# ------------------------------------------------------------
# collapse sample information
# ------------------------------------------------------------

dim(sam.Reinius)
dim(sam.Limbach)
dim(sam.Reynolds.CD4)
dim(sam.Reynolds.mono)
dim(sam.Vento)
dim(sam.Coit)
dim(sam.Zhang)
dim(sam.Schlums)

table(sam.Reinius$id == names(dat.Reinius)[-1])

sams.All = sam.Reinius
studies  = rep("Reinius", nrow(sam.Reinius))

dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

dim(sam.Limbach)
sam.Limbach[1:2,]
table(sam.Limbach$Assay_Name == names(dat.Limbach)[-1])

studies  = c(studies, rep("Limbach", nrow(sam.Limbach)))

sams.All = rbindlist(list(sams.All, sam.Limbach[,.(Assay_Name, cell_type)]))
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

dim(sam.Reynolds.CD4)
sam.Reynolds.CD4[1:2,]
table(sam.Reynolds.CD4$Assay_Name == names(dat.Reynolds.CD4)[-1])

studies  = c(studies, rep("Reynolds", nrow(sam.Reynolds.CD4)))

sams.All = rbindlist(list(sams.All, sam.Reynolds.CD4[,.(Assay_Name, Sample_source_name)]))
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

dim(sam.Reynolds.mono)
sam.Reynolds.mono[1:2,]
table(sam.Reynolds.mono$Assay_Name == names(dat.Reynolds.mono)[-1])

studies  = c(studies, rep("Reynolds", nrow(sam.Reynolds.mono)))

sams.All = rbindlist(list(sams.All, sam.Reynolds.mono[,.(Assay_Name, Sample_source_name)]))
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

dim(sam.Vento)
sam.Vento[1:2,]
table(sam.Vento$Assay_Name == names(dat.Vento)[-1])

studies  = c(studies, rep("Vento", nrow(sam.Vento)))

sams.All = rbindlist(list(sams.All, sam.Vento[,.(Assay_Name, cell_type)]))
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

dim(sam.Coit)
sam.Coit[1:2,]
table(sam.Coit$Assay_Name == names(dat.Coit)[-1])

studies  = c(studies, rep("Coit", nrow(sam.Coit)))

sams.All = rbindlist(list(sams.All, sam.Coit[,.(Assay_Name, cell_type)]))
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

dim(sam.Zhang)
sam.Zhang[1:2,]
table(sam.Zhang$Assay_Name == names(dat.Zhang)[-1])

studies  = c(studies, rep("Zhang", nrow(sam.Zhang)))

sams.All = rbindlist(list(sams.All, sam.Zhang[,.(Assay_Name, cell_type)]))

dim(sam.Schlums)
sam.Schlums[1:2,]
table(sam.Schlums$Assay_Name == names(dat.Schlums)[-1])

studies  = c(studies, rep("Schlums", nrow(sam.Schlums)))

sams.All = rbindlist(list(sams.All, sam.Schlums[,.(Assay_Name, cell_type)]))
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

sams.All[["study"]] = studies
dim(sams.All)
sams.All[1:2,]
table(sams.All$label)

# ------------------------------------------------------------
# re-name cell types
# ------------------------------------------------------------

Tregs = c("act_rTreg", "rTreg")
sams.All$label[which(sams.All$label %in% Tregs)] = "Treg"

DCs = c("activated/mature DCs")
sams.All$label[which(sams.All$label %in% DCs)] = "DC"

Macrophage = c("activated/mature MACs")
sams.All$label[which(sams.All$label %in% Macrophage)] = "Macrophage"

Monocyte = c("CD14+ cell", "CD14+ Monocytes")
Monocyte = c(Monocyte, "peripheral blood isolated CD14+ cells (monocytes)")
sams.All$label[which(sams.All$label %in% Monocyte)] = "Monocyte"

ulabls = unique(sams.All$label)
NK = ulabls[grep("NK", ulabls)]
NK
sams.All$label[which(sams.All$label %in% NK)] = "NK"

B = c("CD19+ B cells")
sams.All$label[which(sams.All$label %in% B)] = "B"

CD4T = ulabls[grep("CD4", ulabls)]
CD4T
sams.All$label[which(sams.All$label %in% CD4T)] = "CD4T"

CD8T = ulabls[grep("CD8", ulabls)]
CD8T
sams.All$label[which(sams.All$label %in% CD8T)] = "CD8T"

Neutrophil = c("Neutrophils", "normal neutrophils")
sams.All$label[which(sams.All$label %in% Neutrophil)] = "Neutrophil"

table(sams.All$label)
table(sams.All$label, sams.All$study)

dim(sams.All)
sams.All[1:5,]

# ------------------------------------------------------------
# combine data
# ------------------------------------------------------------

unique(sams.All$study)

dim(dat.Reinius)
dim(dat.Limbach)
dim(dat.Reynolds.CD4)
dim(dat.Reynolds.mono)
dim(dat.Vento)
dim(dat.Coit)
dim(dat.Zhang)
dim(dat.Schlums)

cpgs = intersect(dat.Reinius$id, dat.Limbach$id)
length(cpgs)

cpgs = intersect(cpgs, dat.Reynolds.CD4$id)
length(cpgs)

cpgs = intersect(cpgs, dat.Reynolds.mono$id)
length(cpgs)

cpgs = intersect(cpgs, dat.Vento$id)
length(cpgs)

cpgs = intersect(cpgs, dat.Coit$id)
length(cpgs)

cpgs = intersect(cpgs, dat.Zhang$id)
length(cpgs)


table(Foxp3 %in% cpgs)
mat1 = match(cpgs, info$ID)
table(info$CHR[mat1])

# since the study of Schlums does not include sex chromsome
# we will not take intersection here. 
# cpgs = intersect(cpgs, dat.Schlums$id)
# length(cpgs)

info1 = info[which(info$ID %in% cpgs),]
dim(info1)
info1[1:2,]

table(Foxp3 %in% info1$ID)

data.All = dat.Reinius[match(info1$ID, dat.Reinius$id),]
data.All = cbind(data.All, dat.Limbach[match(info1$ID, dat.Limbach$id), -1, with=FALSE])
data.All = cbind(data.All, dat.Reynolds.CD4[match(info1$ID, dat.Reynolds.CD4$id), -1, with=FALSE])
data.All = cbind(data.All, dat.Reynolds.mono[match(info1$ID, dat.Reynolds.mono$id), -1, with=FALSE])
data.All = cbind(data.All, dat.Vento[match(info1$ID, dat.Vento$id), -1, with=FALSE])
data.All = cbind(data.All, dat.Coit[match(info1$ID, dat.Coit$id), -1, with=FALSE])
data.All = cbind(data.All, dat.Zhang[match(info1$ID, dat.Zhang$id), -1, with=FALSE])
data.All = cbind(data.All, dat.Schlums[match(info1$ID, dat.Schlums$id), -1, with=FALSE])

dim(data.All)
data.All[1:2,]

table(names(data.All)[-1] == sams.All$id)

# -----------------------------------------------------------------
# select probes that are cell-type-specific
# -----------------------------------------------------------------

dim(sams.All)
sams.All[1:5,]

table(sams.All$label)

ctypes = setdiff(unique(sams.All$label), c("PBMC", "Whole blood"))
wctype = which(sams.All$label %in% ctypes)

table(data.All$id == info1$ID)

sams.ctype = sams.All[wctype,]
data.ctype = data.All[,which(names(data.All) %in% sams.ctype$id), with=FALSE]

dim(data.ctype)
data.ctype[1:2,1:5]

data.ctype = data.matrix(data.ctype)
data.ctype[1:2,1:5]

dim(sams.ctype)
sams.ctype[1:2,]

nna = colSums(is.na(data.ctype))
table(sams.ctype$study[nna > 0])

table(sams.ctype$study)
table(sams.ctype$label)

utypes = unique(sams.ctype$label)
length(utypes)

v   = as.numeric(data.ctype[1,])
lbs = sams.ctype$label

ttest <- function(v, lbs, lb1){
  t1 = t.test(x=v[which(lbs == lb1)], y=v[which(lbs != lb1)])
  t1$p.value
}

pvs = matrix(nrow=nrow(data.ctype), ncol=length(utypes))

for(j in 1:length(utypes)){
  cat(j, date(), "\n")
  lb1 = utypes[j]
  pvs[,j] = apply(data.ctype, 1, ttest, lbs=sams.ctype$label, lb1=lb1)
}

summary(pvs)
colSums(pvs < 0.01)
colSums(pvs < 1e-3)
colSums(pvs < 1e-4)
colSums(pvs < 1e-5)
colSums(pvs < 1e-6)
pcut = 0.05/nrow(info1)
pcut
colSums(pvs < pcut)

pvs[match(Foxp3, info1$ID),]
info1[match(Foxp3, info1$ID),]

pv.min = apply(pvs, 1, min)
summary(pv.min)

table(pv.min < pcut)

dat = data.ctype[which(pv.min < pcut),]
dim(dat)

# ------------------------------------------------------------
# joint PCA
# ------------------------------------------------------------

dat4Pr = dat - rowMeans(dat, na.rm=TRUE)

dat4Pr[is.na(dat4Pr)] = 0
covdat = t(dat4Pr) %*% dat4Pr / nrow(dat4Pr)
dim(covdat)
prdat  = eigen(covdat)

prdat$values[1:20]

PC1 =  prdat$vectors[,1]
PC2 =  prdat$vectors[,2]
PC3 =  prdat$vectors[,3]

sams.ctype$study1 = sams.ctype$study
sams.ctype$study1[which(sams.ctype$study %in% c("Vento", "Zhang"))] = "VZ"
table(sams.ctype$study1)

ev = data.frame(index=1:20, eigenValue=prdat$values[1:20])
p1 <- ggplot(ev, aes(index, eigenValue)) + geom_bar(stat = "identity")
p1

ggsave("../figures/methylation_PCA_eigen_values.pdf", width=9, height=6, 
       units="in")

p2 <- ggplot(sams.ctype, aes(PC1, PC2, shape = factor(study1)))
p2 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_PC1_vs_PC2.pdf", width=9, height=6, 
       units="in")

p3 <- ggplot(sams.ctype, aes(PC1, PC3, shape = factor(study1)))
p3 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_PC1_vs_PC3.pdf", width=9, height=6, 
       units="in")

p4 <- ggplot(sams.ctype, aes(PC2, PC3, shape = factor(study1)))
p4 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_PC2_vs_PC3.pdf", width=9, height=6, 
       units="in")

table(sams.ctype$study1, sams.ctype$label)

# ------------------------------------------------------------
# remove the study of Vento since it looks too different, 
# and remove batch effect at m-value scale
# ------------------------------------------------------------

dim(data.ctype)
dim(sams.ctype)
w2rm = which(sams.ctype$study == "Vento")
length(w2rm)

sams.ctype = sams.ctype[-w2rm,]
data.ctype = data.ctype[,-w2rm]

table(sams.ctype$study, sams.ctype$label)
table(colnames(data.ctype) == sams.ctype$id)

# ------------------------------------------------------------
# remove cell types Eosinophils and Granulocytes
# ------------------------------------------------------------

dim(data.ctype)
dim(sams.ctype)
w2rm = which(sams.ctype$label %in% c("Eosinophils", "Granulocytes"))
length(w2rm)

sams.ctype = sams.ctype[-w2rm,]
data.ctype = data.ctype[,-w2rm]

table(sams.ctype$study, sams.ctype$label)
table(colnames(data.ctype) == sams.ctype$id)

# ------------------------------------------------------------
# remove batch effect at m-value scale
# ------------------------------------------------------------

min(c(data.ctype), na.rm=TRUE)
data.ctype[data.ctype < 1e-4] = 1e-4
min(c(data.ctype), na.rm=TRUE)
c1 = colSums(data.ctype <= 1e-4, na.rm=TRUE)
c1
c2 = rowSums(data.ctype <= 1e-4, na.rm=TRUE)
table(c2 > 0)
table(c2 > 2)
table(c2 > 5)
table(c2 > 10)

table(sams.ctype$study[c1 > 0])
table(sams.ctype$study)
max(c(data.ctype), na.rm=TRUE)

nna = colSums(is.na(data.ctype))
table(sams.ctype$study[nna > 0])

data.ctype.m = log(data.ctype / (1 - data.ctype))

data.ctype.m.Nm = matrix(NA, nrow=nrow(data.ctype.m), ncol=ncol(data.ctype.m))
colnames(data.ctype.m.Nm) = colnames(data.ctype.m)

grps = setdiff(unique(sams.ctype$label), "Treg")
grps

dim(dat.Schlums)
dim(data.ctype)

table(info1$CHR, is.na(data.ctype[,189]))

table(info1$CHR %in% c("X", "Y"), is.na(data.ctype[,189]))

library(preprocessCore)

for(nm1 in grps){
  cat(nm1, date(), "\n")
  
  ww1  = which(sams.ctype$label == nm1 & sams.ctype$study=="Reinius")

  data.ctype.m.Nm[,ww1] = normalize.quantiles(data.ctype.m[,ww1])

  vals = apply(data.ctype.m.Nm[,ww1], 2, sort)

  if(cor(vals[,1], vals[,2]) < 0.99){stop("cor is too small\n")}
  
  vals = rowMeans(vals)

  ww2  = which(sams.ctype$label == nm1)
  if(nm1 == "CD4T"){
    ww2 = c(ww2, which(sams.ctype$label == "Treg"))
  }
  ww2

  for(j in 1:length(ww2)){
    
    if(sams.ctype$study[ww2[j]] == "Schlums"){
      
      xj    = data.ctype.m[,ww2[j]]
      wnaj  = which(!is.na(xj))
      vals1 = vals[wnaj]
      oj    = order(xj[wnaj])
      
      data.ctype.m.Nm[wnaj[oj],ww2[j]]= vals1
    }else{
      oj = order(data.ctype.m[,ww2[j]])
      data.ctype.m.Nm[oj,ww2[j]] = vals
    }
  }
}

nna = colSums(is.na(data.ctype))
table(sams.ctype$study[nna > 0])

nna = colSums(is.na(data.ctype.m.Nm))
table(sams.ctype$study[nna > 0])

dim(data.ctype.m.Nm)
data.ctype.m.Nm[1:2,]

cr1 = cor(data.ctype.m.Nm, data.ctype.m, use="pair")
dr1 = diag(cr1)
summary(dr1)
table(dr1 > 0.97)

diag(cr1) = NA
summary(as.numeric(cr1))
ww3 = which(cr1 > 0.97, TRUE)
dim(ww3)

table(paste(sams.ctype$label[ww3[,1]], sams.ctype$label[ww3[,2]]))

data.ctype = exp(data.ctype.m.Nm)/(1 + exp(data.ctype.m.Nm))

# ------------------------------------------------------------
# select probes
# ------------------------------------------------------------

utypes = unique(sams.ctype$label)
utypes
length(utypes)

ttest <- function(v, lbs, lb1){
  t1 = t.test(x=v[which(lbs == lb1)], y=v[which(lbs != lb1)])
  t1$p.value
}

pvs = matrix(nrow=nrow(data.ctype), ncol=length(utypes))

for(j in 1:length(utypes)){
  cat(j, date(), "\n")
  lb1 = utypes[j]
  pvs[,j] = apply(data.ctype, 1, ttest, lbs=sams.ctype$label, lb1=lb1)
}

summary(pvs)
colSums(pvs < 0.01)
colSums(pvs < 1e-3)
colSums(pvs < 1e-4)
colSums(pvs < 1e-5)
colSums(pvs < 1e-6)
pcut = 0.05/nrow(info1)
pcut
colSums(pvs < pcut)

pvs[match(Foxp3, info1$ID),]
info1[match(Foxp3, info1$ID),]

pv.min = apply(pvs, 1, min)
summary(pv.min)

table(pv.min < pcut)

dat = data.ctype[which(pv.min < pcut),]
dim(dat)
dat[1:2,1:5]

info2 = info1[which(pv.min < pcut),]
dim(info2)
info2[1:2,1:5]

# ------------------------------------------------------------
# joint PCA
# ------------------------------------------------------------

dat4Pr = dat - rowMeans(dat, na.rm=TRUE)

dat4Pr[is.na(dat4Pr)] = 0
covdat = t(dat4Pr) %*% dat4Pr / nrow(dat4Pr)
dim(covdat)
prdat  = eigen(covdat)

prdat$values[1:20]

PC1 =  prdat$vectors[,1]
PC2 =  prdat$vectors[,2]
PC3 =  prdat$vectors[,3]

sams.ctype$study1 = sams.ctype$study
table(sams.ctype$study1)

ev = data.frame(index=1:20, eigenValue=prdat$values[1:20])
p1 <- ggplot(ev, aes(index, eigenValue)) + geom_bar(stat = "identity")
p1

ggsave("../figures/methylation_PCA_after_nm_eigen_values.pdf", width=9, height=6, 
       units="in")

p2 <- ggplot(sams.ctype, aes(PC1, PC2, shape = factor(study1)))
p2 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_after_nm_PC1_vs_PC2.pdf", width=9, height=6, 
       units="in")

p3 <- ggplot(sams.ctype, aes(PC1, PC3, shape = factor(study1)))
p3 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_after_nm_PC1_vs_PC3.pdf", width=9, height=6, 
       units="in")

p4 <- ggplot(sams.ctype, aes(PC2, PC3, shape = factor(study1)))
p4 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_after_nm_PC2_vs_PC3.pdf", width=9, height=6, 
       units="in")

table(sams.ctype$study, sams.ctype$label)

# ------------------------------------------------------------
# regress out study effect
# ------------------------------------------------------------

datNew = dat

for(i in 1:nrow(dat)){
  yi = dat[i,]
  li = lm(yi ~ PC2)
  datNew[i,] = yi - PC2*li$coef[2]
}

dim(datNew)
datNew[1:2,1:5]

datNew[datNew < 0.001] = 0.001
datNew[datNew > 0.999] = 0.999

dat4Pr = datNew - rowMeans(datNew, na.rm=TRUE)

dat4Pr[is.na(dat4Pr)] = 0
covdat = t(dat4Pr) %*% dat4Pr / nrow(dat4Pr)
dim(covdat)
prdat  = eigen(covdat)

prdat$values[1:20]

PC1 =  prdat$vectors[,1]
PC2 =  prdat$vectors[,2]
PC3 =  prdat$vectors[,3]

sams.ctype$study1 = sams.ctype$study
table(sams.ctype$study1)

ev = data.frame(index=1:20, eigenValue=prdat$values[1:20])
p1 <- ggplot(ev, aes(index, eigenValue)) + geom_bar(stat = "identity")
p1

ggsave("../figures/methylation_PCA_after_nm_rmPC2_eigen_values.pdf", 
       width=7.5, height=5, units="in")

p2 <- ggplot(sams.ctype, aes(PC1, PC2, shape = factor(study1)))
p2 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_after_nm_rmPC2_PC1_vs_PC2.pdf", 
       width=7.5, height=5, units="in")

p3 <- ggplot(sams.ctype, aes(PC1, PC3, shape = factor(study1)))
p3 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_after_nm_rmPC2_PC1_vs_PC3.pdf", 
       width=7.5, height=5, units="in")

p4 <- ggplot(sams.ctype, aes(PC2, PC3, shape = factor(study1)))
p4 + geom_point(aes(colour = factor(label)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

ggsave("../figures/methylation_PCA_after_nm_rmPC2_PC2_vs_PC3.pdf", 
       width=7.5, height=5, units="in")

table(sams.ctype$study, sams.ctype$label)

# ------------------------------------------------------------
# write out data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/data")

dat = round(dat, 7)
datNew = round(datNew, 7)

dim(dat)
dat[1:2,1:5]

dim(datNew)
datNew[1:2,1:5]

dim(info2)
info2[1:2,]

dim(sams.ctype)
sams.ctype[1:2,]

sams.ctype = sams.ctype[,1:3, with=FALSE]
dim(sams.ctype)
sams.ctype[1:2,]

write.table(dat, file = "methylation_pure_ct_data.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(datNew, file = "methylation_pure_ct_rmPC2_data.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(sams.ctype, file = "methylation_pure_ct_sample.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(info2, file = "methylation_pure_ct_info.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


quit(save="no")


