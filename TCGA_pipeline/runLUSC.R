
library(data.table)
library(scales)
library(stringr)
library(e1071)
library(ggplot2)
library(MASS)
library(gridExtra)
library(quantreg)
library(quadprog)
library(EMeth)
library(cowplot)

# the raw data of DNA metylation is too large to be kept in gitHub
# here is the local path for DNA methylation data
path.data = "C:/Users/Hanyu/Downloads/Real"
path.work = 'C:/Users/Hanyu/Documents/GitHub/dMeth/TCGA_pipeline'
#source('~/EMeth/EMeth/R/emeth.R')
#source('~/EMeth/EMeth/R/utils.R')
#source('~/EMeth/EMeth/R/cv.emeth.R')
#source('~/EMeth/source/_lib.R')
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
setwd('C:/Users/Hanyu/Downloads/Real')
datM = fread(file = "lusc_methylation.txt", header=TRUE)
sampinfo = fread('lusc_sample_info.txt',header = FALSE)
cpg <- unlist(datM[,1])
colnames(datM) <- gsub('^X',"",colnames(datM))
datM <- datM[,-1]
rownames(datM) = cpg

dir0 = "../TCGA_results/clinical_data/"
setwd('C:/Users/Hanyu/Documents/GitHub/dMeth/TCGA_pipeline/')
tcga_purity = fread(paste0(dir0, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"))
barcode = list()
base_string = '%s-%s-%s'
for(i in 1:length(tcga_purity$array)){
  v1 <- strsplit(tcga_purity$array[i],'-')[[1]][1:3]
  temp <- do.call(sprintf,c(fmt = base_string, as.list(v1)))
  barcode = c(barcode,temp)
}
sampbcr <- sampinfo[,2]
uniq_bcr <- barcode[!(duplicated(barcode)|duplicated(barcode,fromLast = TRUE))]
sampbcr <- sampbcr[!(duplicated(sampbcr)|duplicated(sampbcr,fromLast = TRUE))]
sampinfo <- sampinfo[which(unlist(sampinfo[,2]) %in% unlist(sampbcr)),]
ind <- which(sampinfo$V2 %in% barcode)
mat2 <- unlist(intersect(unlist(sampinfo[,2]),barcode))

patient_id <- sampinfo$V3[match(mat2, sampinfo$V2)] 
datM <- subset(datM, select = patient_id)

purity <- tcga_purity[match(mat2,barcode),]
purity$patient_id <- patient_id
filena <- which(is.na(purity$purity))
purity <- subset(purity, !is.na(purity))
rownames(purity) <- purity$patient_id
if(length(filena) > 0){
datM <- subset(datM, select = -filena)
}
rownames(datM) <- cpg
dim(datM)


# ------------------------------------------------------------
# read in probes to be used
# ------------------------------------------------------------

load(file.path(path.ref, "ref_966probes.RData"))
ls()
length(probe2use)
length(unique(probe2use))

table(probe2use %in% rownames(datM))
table(probe2use %in% rownames(dat))

X  = dat[match(probe2use, rownames(dat)),]
dim(X)
X[1:2,1:5]

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

ys_na      = which(apply(is.na(ys),2,any))
eta_abs_na = which(is.na(purity$purity))

any.na = union(ys_na,eta_abs_na)
any.na

ys = subset(ys,select = -any.na)
purity <- purity[-any.na,]

dim(ys)
ys[1:2,1:5]
table(colnames(ys) == purity$patient_id)
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

#----------------------------------------------------------------------
# Read Estimation from Expression Data, take intersection of the 
# samples with cell type estimation from expression and DNA methylation
#----------------------------------------------------------------------

fnm = '_cibersortx_results/LUSC_composition_cibersortx.txt'
est_expr = fread(fnm)
dim(est_expr)
est_expr[1:2,]

samname  = str_replace(est_expr$Mixture, "^X", "")
length(samname)
samname[1:5]

est_expr = data.matrix(est_expr[,-1])
rownames(est_expr) = samname
dim(est_expr)
est_expr[1:2,]

com_sample = intersect(rownames(est_expr), colnames(ys))
length(com_sample)

est_expr = est_expr[match(com_sample,rownames(est_expr)),]
dim(est_expr)
est_expr[1:2,]

ys     = subset(ys,select = match(com_sample,colnames(ys)))
eta <- purity
eta = eta[match(com_sample,eta$patient_id)]
dim(ys)
ys[1:2,1:4]

table(colnames(ys) == rownames(est_expr))
table(rownames(ys) == rownames(mu))

#----------------------------------------------------------------------
# collapse cell types from expression data into fewer cell types
#----------------------------------------------------------------------

deconv_expr = matrix(NA, nrow = nrow(est_expr),ncol = length(cellTypes))
colnames(deconv_expr) = cellTypes 
rownames(deconv_expr) = rownames(est_expr)
colnames(est_expr)

other = rowSums(est_expr[,c(10,19,20,21)])
deconv_expr[,"B"]    = rowSums(est_expr[,1:3])/0.4
deconv_expr[,"CD4T"] = rowSums(est_expr[,5:8])/0.4
deconv_expr[,"CD8T"] = as.matrix(est_expr[,4])/0.4
deconv_expr[,"Treg"] = as.matrix(est_expr[,9])/0.4
deconv_expr[,"NK"] = rowSums(est_expr[,11:12])/0.42
deconv_expr[,"Monocyte"] = rowSums(est_expr[,13:18])/1.40
deconv_expr[,"Neutrophil"] = as.matrix(est_expr[,22])/0.15
deconv_expr = deconv_expr / rowSums(deconv_expr)

dim(deconv_expr)
deconv_expr[1:2,]

#---------------------------------------------------------------------
# Compare results with different methods and with expression data
#---------------------------------------------------------------------

eta = eta$purity
summary(eta)
eta[which(eta > 0.99)] = 0.99

temp <- rownames(deconv_expr)
deconv_expr <- diag(1-eta) %*% deconv_expr
rownames(deconv_expr) <- temp

mu[mu < 0.05] = 0.05
mu[mu > 0.95] = 0.95

penalty = dim(ys)[1]*(10^seq(-2,1,1)) 

methods = c("EMeth","svr","ls","rls","qp")
rho     = array(data = NA, dim = c(ncol(ys), length(cellTypes), length(methods)),
                dimnames = list(colnames(ys), cellTypes, methods))

alpha   = rep(1/length(cellTypes), length(cellTypes))
simsize = ncol(ys)
ys <- as.matrix(ys)
C = c(0.1,1/sqrt(10),1,sqrt(10),10)

for(j in 1:ncol(ys)){
  if(j %% 10 == 0){ cat(j, date(), "\n") }
  y    = ys[,j]
  X    = as.data.frame(mu)
  Xmat = mu
  
  cv_svr = rep(0,5)
  svrmodel1       = svm(y~., data = X, kernel='linear', cost=0.1, cross=5)
  cv_svr[1]       = mean(svrmodel1$MSE)
  svrmodel2       = svm(y~., data = X, kernel='linear', cost=1/sqrt(10), cross=5)
  cv_svr[2]       = mean(svrmodel2$MSE)
  svrmodel3       = svm(y~., data = X, kernel='linear', cost=1, cross=5)
  cv_svr[3]       = mean(svrmodel2$MSE)
  svrmodel4       = svm(y~., data = X, kernel='linear', cost=sqrt(10), cross=5)
  cv_svr[4]       = mean(svrmodel2$MSE)
  svrmodel5       = svm(y~., data = X, kernel='linear', cost=10, cross=5)
  cv_svr[5]       = mean(svrmodel5$MSE)
  best_svr        = which.min(cv_svr)
  svrmodel        = svm(y~., data = X, kernel='linear', cost=C[best_svr])
  temp            = (t(svrmodel$coefs) %*% svrmodel$SV)
  temp[temp < 0]  = 0
  rho[j,,'svr']   = (1-eta[j])*temp/sum(temp)
  
  temp            = lm(y ~ .-1,data = X)$coefficients
  temp[temp < 0]  = 0
  rho[j,,'ls']    = (1-eta[j])*temp/sum(temp)
  
  temp            = rlm(y ~ .-1,data = X)$coefficients
  temp[temp < 0]  = 0
  rho[j,,'rls']   = (1-eta[j])*temp/sum(temp)
  
  A = rbind(diag(rep(1,length(cellTypes))),rep(-1,length(cellTypes)))
  b = c(rep(0,length(cellTypes)),-1+eta[j])
  D = t(Xmat) %*% Xmat
  d = t(t(Xmat) %*% y)
  rho[j,,'qp']   = (solve.QP(D,d,t(A),b)$solution)
}

print('EMeth')
hundrediter = cv.emeth(ys, eta, mu, aber = TRUE, V='c', init = 'default',
                       family = 'normal', nu = penalty, folds = 5, 
                       maxiter = 50, verbose = TRUE)
rho[,,'EMeth'] = hundrediter[[1]]$rho

#---------------------------------------------------------------------
# save the results
#---------------------------------------------------------------------

dim(rho)
rho[1,,]
dimnames(rho)

dim(deconv_expr)
deconv_expr[1:2,]

rho_LUSC = rho
deconv_expr_LUSC = deconv_expr
save(rho_LUSC, file = '../TCGA_results/deconv_methy_LUSC.RData')
save(deconv_expr_LUSC, file = '../TCGA_results/deconv_expr_LUSC.RData')

#---------------------------------------------------------------------
# generate plots
#---------------------------------------------------------------------

setwd("_figures_LUSC")

utypes = intersect(cellTypes,colnames(deconv_expr))
utypes

cormat <- matrix(NA,nrow = length(utypes), ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- utypes

err <- matrix(NA,nrow = length(utypes), ncol = length(methods))
colnames(err) <- methods
rownames(err) <- utypes

rss <- matrix(NA,nrow = length(utypes), ncol = length(methods))
colnames(rss) <- methods
rownames(rss) <- utypes

for(i in 1:length(utypes)){
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){
    cor(rho[,utypes[i],methods[j]],deconv_expr[,utypes[i]])
  })
  err[i,] <- sapply(1:length(methods), FUN = function(j){
    sqrt(mean((rho[,utypes[i],methods[j]] - deconv_expr[,utypes[i]])^2))
  }) 
  rss[i,] <- sapply(1:length(methods), FUN = function(j){
    temp <- lm(deconv_expr[,utypes[i]]~rho[,utypes[i],methods[j]])
    return(sum(temp$residuals^2))
  })
}

for(i in 1:length(utypes)){
  pdf(sprintf('%s_express_methy_discard_high_purity.pdf',utypes[i]))
  plist = list()
  plist <- lapply(1:length(methods), FUN = function(j){
    tempdata = cbind(rho[,utypes[i],methods[j]],deconv_expr[,utypes[i]],eta)
    colnames(tempdata) <- c("methylation","expression","eta")
    newplot <- ggplot(data = as.data.frame(tempdata), 
                      aes(x=methylation,y=expression,color=eta)) + 
      xlim(0,0.3) + ylim(0,0.3) + geom_point() + 
      geom_abline(intercept = 0,slope = 1) + ggtitle(sprintf('%s on %s',methods[j],utypes[i]))
  })
  grid.arrange(grobs = plist,ncol=2)
  save(plist, file = sprintf('plist_%s_%s.RData',utypes[i],'LUSC'))
  dev.off()
}

pdf('correlation.pdf')
tempdata <- melt(as.data.table(cormat))
colnames(tempdata) <- c('Methods','Correlation')
tempdata$cellType = rep(utypes,5)
p1 <- ggplot(tempdata,aes(x=Methods,y=Correlation)) + geom_boxplot() +
  geom_point(size = 5,aes(colour = cellType)) + theme_cowplot() + ggtitle('Correlation for LUSC')
print(p1)
save(p1,file=sprintf('Cor_%s.RData','LUSC'))
dev.off()


pdf('correlation.pdf')
tempdata <- melt(as.data.table(err))
colnames(tempdata) <- c('Methods','RMSE')
tempdata$cellType = rep(utypes,5)
p2 <- ggplot(tempdata,aes(x=Methods,y=RMSE)) + geom_boxplot() +
  geom_point(size = 5,aes(colour = cellType)) + theme_cowplot() + ggtitle('RMSE for LUSC')
print(p2)
save(p2,file=sprintf('err_%s.RData','LUSC'))
dev.off()


print(cormat)
print(err)
print(rss)

OneMinusCorr <- 1-matrix(cormat,ncol = 1, byrow = FALSE)
RMSE <- matrix(err,ncol = 1, byrow = FALSE )
CellType <- rep(cellTypes,length(methods))
Methods <- rep(methods,each = length(cellTypes))
res <- cbind.data.frame(OneMinusCorr,RMSE,CellType,Methods)

pdf('Comparison.pdf')
complot<- ggplot(res, aes(x=OneMinusCorr,y=RMSE, color =  Methods)) + 
  ggtitle('LUSC') + geom_point() + 
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks('log10',function(x) 10^x),
                     labels = trans_format('log10',math_format(10^.x))) + 
  geom_text(label = res[,3])+xlim(0.1,1.05)
print(complot)
dev.off()
cormat_LUSC <- cormat
err_LUSC <- err
save(cormat_LUSC, file='cormat_LUSC.RData')
save(err_LUSC, file = 'err_LUSC.RData')

sessionInfo()
gc()

quit(save = 'no')



