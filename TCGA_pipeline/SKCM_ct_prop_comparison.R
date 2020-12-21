
library(data.table)
library(stringr)
library(ggplot2)
library(MASS)
library(gridExtra)
library(cowplot)
library(ggcorrplot)

# ------------------------------------------------------------
# read in cell type proportion estimates by EMeth
# ------------------------------------------------------------

load("TCGA_results/deconv_methy_SKCM.RData")
ls()
dim(rho_SKCM)
rho_SKCM[1,,]

# ------------------------------------------------------------
# read in cell type proportion estimates by CYBERSORT
# cb1 was generatd using count as input
# cb2 was generated using TPM as input
# cb3 was generated using TPM and a new reference from SF2018
# ------------------------------------------------------------

rho_cb1 = fread("_cibersortx_results/SKCM_composition_cibersortx.txt")
dim(rho_cb1)
rho_cb1[1:2,1:5]

rho_cb2 = fread("_cibersortx_results/CIBERSORTx_LM22_Adjusted.txt")
dim(rho_cb2)
rho_cb2[1:2,1:5]

rho_cb3 = fread("_cibersortx_results/CIBERSORTx_SF2018_Adjusted.txt")
dim(rho_cb3)
rho_cb3[1:2,1:5]

names(rho_cb1) = gsub(" ", ".", names(rho_cb1), fixed=TRUE)
table(names(rho_cb1) == names(rho_cb2))

w1 = which(names(rho_cb1) != names(rho_cb2))
names(rho_cb1)[w1]
names(rho_cb2)[w1]
names(rho_cb1)[w1] = names(rho_cb2)[w1] = "T.cells.regulatory"

rho_cb2$sample  = rho_cb2$Mixture
rho_cb3$sample  = rho_cb3$Mixture
rho_cb2$Mixture = substr(rho_cb2$Mixture, 9, 12)
rho_cb3$Mixture = substr(rho_cb3$Mixture, 9, 12)

table(rho_cb1$Mixture %in% rho_cb2$Mixture)
table(rho_cb1$Mixture %in% rho_cb3$Mixture)

rho_cb2 = rho_cb2[match(rho_cb1$Mixture, rho_cb2$Mixture),]
dim(rho_cb2)
rho_cb2[1:2,]

rho_cb3 = rho_cb3[match(rho_cb1$Mixture, rho_cb3$Mixture),]
dim(rho_cb3)
rho_cb3[1:2,]

cr12 = cor(rho_cb1[,2:23], rho_cb2[,2:23])
dim(cr12)
cr12[1:3,1:3]

cr13 = cor(rho_cb1[,2:23], rho_cb3[,2:12])
dim(cr13)
cr13[1:3,1:3]

cr23 = cor(rho_cb2[,2:23], rho_cb3[,2:12])
dim(cr23)
cr23[1:3,1:3]

ggcorrplot(cr12)
ggcorrplot(cr13)
ggcorrplot(cr23)

# ------------------------------------------------------------
# read in tumor purity information
# ------------------------------------------------------------

dir0 = "TCGA_results/clinical_data/"
tcga_purity = fread(paste0(dir0, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"))

dim(tcga_purity)
tcga_purity[1:2,]

rho_cb2$array = substr(rho_cb2$sample, 1, 15)
rho_cb2$array = gsub(".", "-", rho_cb2$array, fixed=TRUE)
dim(rho_cb2)
rho_cb2[1:2,21:28]

table(rho_cb2$array %in% tcga_purity$array)

mat1 = match(rho_cb2$array, tcga_purity$array)
wnNA = which(!is.na(mat1))
rho_cb2$purity = rep(NA, nrow(rho_cb2))
rho_cb2$purity[wnNA] = tcga_purity$purity[mat1[wnNA]]
table(rho_cb2$array[wnNA] == tcga_purity$array[mat1[wnNA]])

#----------------------------------------------------------------------
# collapse cell types from expression data into fewer cell types
#----------------------------------------------------------------------

ct_list = list()
ct_list[["B"]] = c("B.cells.naive", "B.cells.memory", "Plasma.cells")
ct_list[["CD4T"]] = c("T.cells.CD4.naive", "T.cells.CD4.memory.resting", 
                      "T.cells.CD4.memory.activated", 
                      "T.cells.follicular.helper")
ct_list[["CD8T"]] = c("T.cells.CD8")
ct_list[["Treg"]] = c("T.cells.regulatory")
ct_list[["NK"]]   = c("NK.cells.resting", "NK.cells.activated")
ct_list[["Monocyte"]]   = c("Monocytes", "Macrophages.M0", "Macrophages.M1", 
                            "Macrophages.M2", "Dendritic.cells.resting",
                            "Dendritic.cells.activated")
ct_list[["Neutrophil"]] = c("Neutrophils")

cell_size_factors = c(0.4, 0.4, 0.4, 0.4, 0.42, 1.40, 0.15)
names(cell_size_factors) = c("B", "CD4T", "CD8T", "Treg", "NK", 
                             "Monocyte", "Neutrophil")


rho_expr_lm22 = matrix(NA, nrow = nrow(rho_cb2), ncol = length(ct_list))
colnames(rho_expr_lm22) = names(ct_list) 
rownames(rho_expr_lm22) = rho_cb2$Mixture

for(ct1 in names(ct_list)){
  cts = ct_list[[ct1]]
  rho_expr_lm22[,ct1] = rowSums(as.matrix(rho_cb2[,..cts]))
}

rho_expr_lm22 = rho_expr_lm22 / rowSums(rho_expr_lm22)

dim(rho_expr_lm22)
rho_expr_lm22[1:2,]

#----------------------------------------------------------------------
# collapse cell types from expression data into fewer cell types
#----------------------------------------------------------------------

ct_list = list()
ct_list[["B"]] = c("B_cells", "Plasma_cells")
ct_list[["CD4T"]] = c("CD4T_memory")
ct_list[["CD8T"]] = c("CD8T_B", "CD8T_G")
ct_list[["Treg"]] = c("Tregs")
ct_list[["NK"]]   = c("NK1", "NK2", "NK3")
ct_list[["Monocyte"]]   = c("Monocytes_Macrophages", "Dendritic_cells")
cell_size_factors = c(0.4, 0.4, 0.4, 0.4, 0.42, 1.40)
names(cell_size_factors) = c("B", "CD4T", "CD8T", "Treg", "NK", "Monocyte")

rho_expr_sf11 = matrix(NA, nrow = nrow(rho_cb3), ncol = length(ct_list))
colnames(rho_expr_sf11) = names(ct_list) 
rownames(rho_expr_sf11) = rho_cb3$Mixture

for(ct1 in names(ct_list)){
  cts = ct_list[[ct1]]
  rho_expr_sf11[,ct1] = rowSums(as.matrix(rho_cb3[,..cts]))
}

rho_expr_sf11 = rho_expr_sf11 / rowSums(rho_expr_sf11)

dim(rho_expr_sf11)
rho_expr_sf11[1:2,]

#----------------------------------------------------------------------
# correct for tumor purity
#----------------------------------------------------------------------

table(rownames(rho_expr_lm22) == rownames(rho_expr_sf11))

table(dimnames(rho_SKCM)[[1]] %in% rownames(rho_expr_lm22))

match_meth2expr = match(dimnames(rho_SKCM)[[1]], rownames(rho_expr_lm22))
rho_expr_lm22 = rho_expr_lm22[match_meth2expr,]
rho_expr_sf11 = rho_expr_sf11[match_meth2expr,]

dim(rho_expr_lm22)
dim(rho_expr_sf11)

eta = rho_cb2$purity[match_meth2expr]
summary(eta)
eta[which(eta > 0.99)] = 0.99

rho_expr_lm22 = diag(1-eta) %*% rho_expr_lm22
rho_expr_sf11 = diag(1-eta) %*% rho_expr_sf11

rownames(rho_expr_lm22) = dimnames(rho_SKCM)[[1]]
rownames(rho_expr_sf11) = dimnames(rho_SKCM)[[1]]

# ------------------------------------------------------------
# compare cell type proportion
# ------------------------------------------------------------

methods = c("EMeth","svr","ls","rls","qp")
cellTypes = dimnames(rho_SKCM)[[2]]

utypes = intersect(cellTypes, colnames(rho_expr_lm22))
utypes

cormat <- matrix(NA,nrow = length(utypes), ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- utypes

err = cormat

for(i in 1:length(utypes)){
  rho_e_i = rho_expr_lm22[,utypes[i]]
  
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){
    cor(rho_SKCM[,utypes[i],methods[j]], rho_e_i, use="pair")
  })
  err[i,] <- sapply(1:length(methods), FUN = function(j){
    sqrt(mean((rho_SKCM[,utypes[i],methods[j]] - rho_e_i)^2, na.rm=TRUE))
  }) 
}

cormat
err

# ------------------------------------------------------------
# compare cell type proportion
# ------------------------------------------------------------

utypes = intersect(cellTypes, colnames(rho_expr_sf11))
utypes

cormat <- matrix(NA,nrow = length(utypes), ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- utypes

err = cormat

for(i in 1:length(utypes)){
  rho_e_i = rho_expr_sf11[,utypes[i]]
  
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){
    cor(rho_SKCM[,utypes[i],methods[j]], rho_e_i, use="pair")
  })
  err[i,] <- sapply(1:length(methods), FUN = function(j){
    sqrt(mean((rho_SKCM[,utypes[i],methods[j]] - rho_e_i)^2, na.rm=TRUE))
  }) 
}

cormat
err

Treg_e = rho_expr_sf11[,"Treg"]
Treg_m = rho_SKCM[, "Treg", "EMeth"]

sessionInfo()
gc()

quit(save = 'no')



