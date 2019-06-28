
# ------------------------------------------------------------
# binary classification for gene expression data
#
# select probes with cell-type p-value < pcut1 and R2 > R2cut1
# ------------------------------------------------------------

bclassE <- function(level1, sam, dat, path, R2cut1=0.2){
  
  sam2use = which(sam$label %in% unlist(level1))
  nms     = names(level1)
  
  ctype1  = rep("A", nrow(sam))
  ctype1[which(sam$label %in% level1[[nms[1]]])] = "B"
  ctype1[which(sam$label %in% level1[[nms[2]]])] = "C"

  pvs = R2s = rep(NA, nrow(dat))
  
  w2kp = which(ctype1 %in% c("B", "C"))
  dat1 = dat[,w2kp]
  
  ctype2 = ctype1[w2kp]
  nn = length(w2kp)
  wB = which(ctype2 == "B")
  wC = which(ctype2 == "C")
  
  nB = length(which(ctype2 == "B"))
  nC = length(which(ctype2 == "C"))
  
  vAll = apply(dat1, 1, var, na.rm=TRUE)*(nn - 1)
  varB = apply(dat1[,wB], 1, var, na.rm=TRUE)*(nB - 1)
  varC = apply(dat1[,wC], 1, var, na.rm=TRUE)*(nC - 1)
  R2s  = 1 - (varB + varC)/vAll
  
  w2use = which(R2s > R2cut1)
  
  db1 = data.frame(id=rownames(dat)[w2use], R2=R2s[w2use],
                    stringsAsFactors=FALSE)
  
  db1 = db1[order(db1$R2, decreasing = TRUE),]
  
  if(nrow(db1) > 1000){
    db1 = db1[1:1000,]
  }
  
  db1
}

# ------------------------------------------------------------
# calculate the overlap of two densities
# ------------------------------------------------------------

iOverlap <- function(d0, d1){
  
  xx0 = d0$x
  xx1 = d1$x
  
  if(any(xx0 != xx1)){
    stop("some values of xx0 and xx1 are different\n")
  }
  
  ymax = pmax(d0$y, d1$y)
  
  xl = xx0[-length(xx0)]
  xr = xx0[-1]
  
  yl = ymax[-length(xx0)]
  yr = ymax[-1]
  
  2 - sum((xr - xl)*(yr + yl)/2)
}

# ------------------------------------------------------------
# plot the density of the probes to be used
#
# which has p-value for Shapiro test > pcut, and
# area of overlpa < ocut
# ------------------------------------------------------------

plot1 <- function(level1, dat, sam, rids, pcut=0.1, ocut=0.5, xlab="expression")
{
  nms     = names(level1)
  sam2use = which(sam$label %in% unlist(level1))
  
  ctype1  = rep(2, nrow(sam))
  ctype1[which(sam$label %in% level1[[nms[1]]])] = 0
  ctype1[which(sam$label %in% level1[[nms[2]]])] = 1
  
  table(ctype1)
  
  pvs  = matrix(NA, nrow=length(rids), ncol=10)
  cnm0 = paste0(c("normal_test_", "mean_", "sd_"), nms[1])
  cnm1 = paste0(c("normal_test_", "mean_", "sd_"), nms[2])
  cnm2 = paste0(c("normal_test_", "mean_", "sd_"), "others")
  colnames(pvs) = c(cnm0, cnm1, cnm2, "overlap")
  
  ww0 = which(ctype1==0)
  ww1 = which(ctype1==1)
  ww2 = which(ctype1==2)

  pp0 = length(ww0)/nrow(sam)
  pp1 = length(ww1)/nrow(sam)
  pp2 = length(ww2)/nrow(sam)

  for(i in 1:length(rids)){
    ridi = rids[i]
    wwi  = match(ridi, rownames(dat))
    yi   = dat[wwi,]
    
    min1 = min(yi, na.rm=TRUE)
    max1 = max(yi, na.rm=TRUE)
    
    yi0 = na.omit(yi[ww0])
    yi1 = na.omit(yi[ww1])
    yi2 = na.omit(yi[ww2])

    if(length(unique(yi0)) >= 3){
      d0 = density(yi0, from=min1, to=max1, n=1024)
      t0 = shapiro.test(yi0)
      pvs[i,1] = t0$p.value
      pvs[i,2] = mean(yi0)
      pvs[i,3] = sd(yi0)
    }
    
    if(length(unique(yi1)) >= 3){
      d1 = density(yi1, from=min1, to=max1, n=1024)
      t1 = shapiro.test(yi1)
      pvs[i,4] = t1$p.value
      pvs[i,5] = mean(yi1)
      pvs[i,6] = sd(yi1)
    }

    if(length(unique(yi2)) >= 3){
      d2 = density(yi2, from=min1, to=max1, n=1024)
      t2 = shapiro.test(yi2)
      pvs[i,7] = t2$p.value
      pvs[i,8] = mean(yi2)
      pvs[i,9] = sd(yi2)
    }

    if(length(unique(yi0)) >= 3 && length(unique(yi1)) >= 3){
      pvs[i,10] = iOverlap(d0, d1)
    }
    
    if(any(is.na(pvs[i,c(1,4)]))){ next }
    
    if(all(pvs[i,c(1,4,7)] > pcut, na.rm=TRUE) & pvs[i,10] < ocut){
      
      if(is.na(pvs[i,7])){
        dby = d0$y*pp0 + d1$y*pp1
        ymax = max(c(d0$y, d1$y), na.rm=TRUE)*1.25

      }else{
        dby  = d0$y*pp0 + d1$y*pp1 + d2$y*pp2
        ymax = max(c(d0$y, d1$y, d2$y), na.rm=TRUE)*1.25

      }
      
      
      if(is.na(pvs[i,7])){
        mmi  = sprintf("%s, p.norm=(%.1e, %.1e)", ridi, pvs[i,1], pvs[i,4])
      }else{
        mmi  = sprintf("%s, p.norm=(%.1e, %.1e, %.1e)", ridi, pvs[i,1], pvs[i,4], pvs[i,7])
      }
      
      plot(d0$x, dby, type="n", lwd=2, col="grey", main="", bty="n",
      xlab=xlab, ylab="density", xlim=c(0,1), ylim=c(0,ymax))
      mtext(mmi, line=0.5)
      
      lines(d0$x, d0$y, col="red")
      lines(d1$x, d1$y, col="skyblue")
      
      if(is.na(pvs[i,7])){
        legend("topright", c(names(level1)), lty=c(1,1),
          col=c("red", "skyblue"), bty="n")
      }else{
        legend("topright", c(names(level1), "others"), lty=c(1,1,1),
        col=c("red", "skyblue", "darkgrey"), bty="n")
        lines(d2$x, d2$y, col="darkgrey")
      }

    }
  }
  
  pvs
  
}

# ------------------------------------------------------------
# select probes
# ------------------------------------------------------------

probeSelect <- function(level1, sam, dat, path, dataType="expression"){
  
  if(dataType == "expression"){
    xlab = "gene expression"
  }else if(dataType == "methylation"){
    xlab = "beta-value"
  }
  
  db1  = bclassE(level1, sam, dat, path)

  if(nrow(db1) == 0){ return(NULL) }

  dim(db1)
  db1[1:2,]
  
  nms = names(level1)
  fn  = sprintf("%s/ex_%s_vs_%s.pdf", path, nms[1], nms[2])
  
  pdf(fn, width=5, height=4)
  pvs = plot1(level1, dat, sam, db1$id, xlab=xlab)
  dev.off()
  
  db1 = cbind(db1, pvs)
  dim(db1)
  db1[1:2,]
  
  pvs = pvs[,c(1,4,7), drop=FALSE]
  
  if(setequal(sam$label, unlist(level1))){
    
    if(all(is.na(pvs[,1])) && all(is.na(pvs[,2])) ){
      w2kp = 1:nrow(db1)
    }else {
      if(all(is.na(pvs[,1])) ){
        w2kp = which(pvs[,2] > 0.05)
      }else if(all(is.na(pvs[,2])) ){
        w2kp = which(pvs[,1] > 0.05)
      }else{
        w2kp = which(pvs[,1] > 0.05 & pvs[,2] > 0.05)
      }
    }

  }else{

    if(any(is.na(pvs[,3]))){
      stop("hm, I do not expect normality test to be NA for group 'others'\n")
    }

    if(all(is.na(pvs[,1])) && all(is.na(pvs[,2])) ){
      w2kp = 1:nrow(db1)
    }else {
      if(all(is.na(pvs[,1])) ){
        w2kp = which(pvs[,2] > 0.05 & pvs[,3] > 0.05)
      }else if(all(is.na(pvs[,2])) ){
        w2kp = which(pvs[,1] > 0.05 & pvs[,3] > 0.05)
      }else{
        w2kp = which(pvs[,1] > 0.05 & pvs[,2] > 0.05 & pvs[,3] > 0.05)
      }
    }
  }
  
  if(length(w2kp) == 0){ return(NULL) }
  
  db1  = db1[w2kp,]
  dim(db1)
  db1[1:2,]
  
  #   # ANOVA within each sub-category of cell types
  #
  #   pvs = matrix(NA, nrow=nrow(db1), ncol=2)
  #
  #   wRow = match(db1$id, rownames(dat))
  #
  #   if(any(is.na(wRow))){
  #     stop("some row ids are missing\n")
  #   }
  #
  #   nms = names(level1)
  #
  #   wCol1 = which(sam$label %in% level1[[nms[1]]])
  #   wCol2 = which(sam$label %in% level1[[nms[2]]])
  #
  #   dat1  = dat[wRow, wCol1]
  #   dat2  = dat[wRow, wCol2]
  #
  #   xi1   = as.factor(sam$label[wCol1])
  #   xi2   = as.factor(sam$label[wCol2])
  #
  #   for(i in 1:nrow(db1)){
  #     if(length(level1[[nms[1]]]) > 1){
  #       yi1 = dat1[i,]
  #       li1 = lm(yi1 ~ xi1)
  #       ai1 = anova(li1)
  #       pvs[i,1] = ai1$Pr[1]
  #     }
  #
  #     if(length(level1[[nms[2]]]) > 1){
  #       yi2 = dat2[i,]
  #       li2 = lm(yi2 ~ xi2)
  #       ai2 = anova(li2)
  #       pvs[i,2] = ai2$Pr[1]
  #     }
  #   }
  #
  #   pmin1 = pmin(pvs[,1], pvs[,2], na.rm = TRUE)
  #
  #   if(! all(is.na(pmin1))){
  #
  #     if(any(is.na(pmin1))){
  #       stop("I expect anova p-value to be missing for all or none\n")
  #     }
  #
  #     table(pvs[,1] > 0.01, pvs[,2] > 0.01)
  #     table(pmin1 > 0.001)
  #
  #     table(pmin1 >= 1e-3, db1$pv < 0.01*pmin1)
  # 
  #     db1 = db1[which(pmin1 >= 1e-3 & db1$pv < 0.01*pmin1),]
  #   }
  
  db1
}

# ------------------------------------------------------------
# normal MLE for two cell types
# ------------------------------------------------------------

getCT <- function(level1, sam){
  cts = rep(NA, nrow(sam))
  nms = names(level1)
  
  ww1 = which(sam$label %in% level1[[nms[1]]])
  ww2 = which(sam$label %in% level1[[nms[2]]])
  ww3 = setdiff(1:nrow(sam), c(ww1, ww2))
  
  cts[ww1] = nms[1]
  cts[ww2] = nms[2]
  
  if(length(ww3) > 0)
  cts[ww3] = "others"
  
  cts
}

normal.mixture.binary <- function(rho1, y, X, level1, sam){
  
  cellType = getCT(level1, sam)
  table(cellType)

  if(min(as.numeric(table(cellType))) < 2){
    stop("at least two samples are required each cell type\n")
  }
  
  fun1 <- function(v1, t1){
    tapply(v1, t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(v1, t1, sd)
  }

  muX  = t(apply(X, 1, fun1, t1=cellType))
  sdX  = t(apply(X, 1, fun2, t1=cellType))

  if(ncol(muX) != 2){
    stop("I expect binary mixture\n")
  }
  
  mat1 = match(names(level1), colnames(muX))
  muX  = muX[,mat1]
  sdX  = sdX[,mat1]
  
  muY  = muX %*% c(rho1, 1 -rho1)
  vrY  = sdX^2 %*% c(rho1^2, (1 -rho1)^2)
  
  # sdY = sqrt(sum(res^2)/length(y))
  # -log(sdY)
  
  res = y - muY
  # z   = res/sqrt(vrY)
  # sum(dnorm(z, log=TRUE))
  -sum(res^2)
}

# ------------------------------------------------------------
# nnls for multiple cell types
# ------------------------------------------------------------

wls <- function(y, X, cellTypes, sam, total=1.0){
  
  wSam = which(sam$label %in% cellTypes)
  sam1 = sam[wSam,]
  X1   = X[,wSam]
  
  wnNA = which(!is.na(y))
  wnNA = intersect(wnNA, which(rowSums(is.na(X1)) == 0))
  y1   = y[wnNA]
  X1   = X1[wnNA,]
  
  fun1 <- function(v1, t1){
    tapply(v1, t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(v1, t1, sd)
  }
  
  muX  = t(apply(X1, 1, fun1, t1=sam1$label))
  sdX  = t(apply(X1, 1, fun2, t1=sam1$label))
  
  muX  = muX[,match(cellTypes, colnames(muX))]
  sdX  = sdX[,match(cellTypes, colnames(sdX))]

  muX1 = cbind(rep(1, nrow(muX)), muX)
  
  n0 = nnls(muX1, y1)
  b0 = coef(n0)[-1]
  b0 = total*b0/sum(b0)
  b1 = b0
  
  for(kk in 1:100){
    sdY  = drop(sqrt(sdX^2 %*% b1^2))
    y2   = y1/sdY
    muX2 = muX1/sdY
    n2 = nnls(muX2, y2)
    b2 = coef(n2)[-1]
    b2 = total*b2/sum(b2)
    
    if(max(abs(b1 - b2)) < 1e-5){
      break
    }
    b1 = b2
  }
  
  list(cellType=cellTypes, b0=b0, b2=b2, nIt=kk)
}

# ------------------------------------------------------------
# EM algorithm for mixture of a normal distribution
# plus a triangluar distribution
# ------------------------------------------------------------

mixEM <- function(r){
  
  if(any(abs(r) > 0.99)){
    warning("values with abs > .99 are truncated\n")
    r[which(r > 0.99)]  = 0.99
    r[which(r < -0.99)] = -0.99
  }
  
  sd1  = sd(r)
  nn1  = length(r)
  den1 = dnorm(r, mean=0, sd=sd1)
  den2 = 1-abs(r)
  alpha = 0.8
  
  for(k in 1:100){
    sd0  = sd1

    pp1  = alpha*den1/(alpha*den1 + (1-alpha)*den2)
    sd1  = sqrt(sum(r*r*pp1)/(nn1 - 1))
    den1 = dnorm(r, mean=0, sd=sd1)
    alpha = sum(pp1)/nn1

    if(abs(sd0 - sd1) < 1e-5){break}
    
    if(sd1 < 0.01){break}
  }
  
  list(p1=drop(pp1), alpha=alpha, sd1=sd1)
}

# ------------------------------------------------------------
# nnls for multiple cell types with errors
# ------------------------------------------------------------

wlsEM <- function(y, X, cellTypes, sam, total=1.0){
  
  wSam = which(sam$label %in% cellTypes)
  sam1 = sam[wSam,]
  X1   = X[,wSam]
  
  wnNA = which(!is.na(y))
  wnNA = intersect(wnNA, which(rowSums(is.na(X1)) == 0))
  y1   = y[wnNA]
  X1   = X1[wnNA,]
  
  fun1 <- function(v1, t1){
    tapply(v1, t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(v1, t1, sd)
  }
  
  muX  = t(apply(X1, 1, fun1, t1=sam1$label))
  sdX  = t(apply(X1, 1, fun2, t1=sam1$label))
  
  muX  = muX[,match(cellTypes, colnames(muX))]
  sdX  = sdX[,match(cellTypes, colnames(sdX))]
  
  muX1 = cbind(rep(1, nrow(muX)), muX)
  
  n0 = nnls(muX1, y1)
  e1 = mixEM(n0$resid)

  b0 = coef(n0)[-1]
  b0 = total*b0/sum(b0)
  b1 = b0
  
  for(kk in 1:100){
    sdY  = drop(sqrt(sdX^2 %*% b1^2))
    wts  = e1$p1/sdY
    y2   = y1*wts
    muX2 = muX1*wts
    n2 = nnls(muX2, y2)
    b2 = coef(n2)[-1]
    b2 = total*b2/sum(b2)
    
    if(max(abs(b1 - b2)) < 1e-5){
      break
    }

    b1 = b2
    
    r1 = n2$resid/wts
    e1 = mixEM(n0$resid)
  }
  
  list(cellType=cellTypes, b0=b0, b2=b2, nIt=kk, nobs=sum(e1$p1))
}

# ------------------------------------------------------------
# elastic net for multiple cell types
# ------------------------------------------------------------

enet <- function(y, X, cellTypes, sam, alpha, total=1.0){
  
  wSam = which(sam$label %in% cellTypes)
  sam1 = sam[wSam,]
  X1   = X[,wSam]
  
  wnNA = which(!is.na(y))
  wnNA = intersect(wnNA, which(rowSums(is.na(X1)) == 0))
  y1   = y[wnNA]
  X1   = X1[wnNA,]
  
  fun1 <- function(v1, t1){
    tapply(v1, t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(v1, t1, sd)
  }
  
  muX  = t(apply(X1, 1, fun1, t1=sam1$label))
  sdX  = t(apply(X1, 1, fun2, t1=sam1$label))
  
  muX  = muX[,match(cellTypes, colnames(muX))]
  sdX  = sdX[,match(cellTypes, colnames(sdX))]
  
  cg1  = cv.glmnet(muX, y1, alpha=alpha)
  b1   = coef(cg1, s = "lambda.min")[-1]
  b1[which(b1 < 0)] = 0.01
  b1 = total*b1/sum(b1)
  b0 = b1
  
  for(kk in 1:100){
    sdY  = drop(sqrt(sdX^2 %*% b1^2))
    cg2  = cv.glmnet(muX, y1, weights=1/sdY, alpha=alpha)
    b2   = coef(cg2, s = "lambda.min")[-1]
    b2[which(b2 < 0)] = 0.01
    b2   = total*b2/sum(b2)
    
    if(max(abs(b1 - b2)) < 1e-5){
      break
    }

    b1   = b2
  }
  
  list(cellType=cellTypes, b0=b0, b2=b2, nIt=kk)
  
}

# ------------------------------------------------------------
# logNormal MLE for two cell types
# ------------------------------------------------------------

log.normal.mixture.binary <- function(para, y, x1, x2){
  
  rho1   = para[1]
  delta1 = para[2]
  
  n1  = ncol(x1)
  n2  = ncol(x2)
  
  if(n1 < 2 || n2 < 2){
    stop("at least two samples are required for x1 and x2\n")
  }
  
  logX1 = log(x1)
  logX2 = log(x2)
  
  mu1 = apply(logX1, 1, mean) + log(rho1*delta1)
  mu2 = apply(logX2, 1, mean) + log((1 - rho1)*delta1)
  
  sd1 = apply(logX1, 1, sd)
  sd2 = apply(logX2, 1, sd)
  
  muX     = cbind(mu1, mu2)
  sigma2X = cbind(sd1*sd1, sd2*sd2)
  
  ems     = exp(muX + sigma2X/2)
  sigma2Y = log(rowSums(ems*ems*(exp(sigma2X) -1))/(rowSums(ems))^2 + 1)
  muY     = log(rowSums(ems)) - sigma2Y/2

  logY = log(y)
  zY   = (logY - muY)/sqrt(sigma2Y)
  
  sum(dnorm(zY, log=TRUE))
}


log.normal.mixture.knownScale <- function(rho1, delta1, y, x1, x2){
  
  n1  = ncol(x1)
  n2  = ncol(x2)
  
  if(n1 < 2 || n2 < 2){
    stop("at least two samples are required for x1 and x2\n")
  }
  
  logX1 = log(x1)
  logX2 = log(x2)
  
  mu1 = apply(logX1, 1, mean)
  mu2 = apply(logX2, 1, mean)


  mu1 = apply(logX1, 1, mean) + log(rho1*delta1)
  mu2 = apply(logX2, 1, mean) + log((1 - rho1)*delta1)
  
  sd1 = apply(logX1, 1, sd)*sqrt((n1 - 1)/n1)
  sd2 = apply(logX2, 1, sd)*sqrt((n2 - 1)/n2)
  
  muX     = cbind(mu1, mu2)
  sigma2X = cbind(sd1*sd1, sd2*sd2)
  
  ems     = exp(muX + sigma2X/2)
  sigma2Y = log(rowSums(ems*ems*(exp(sigma2X) -1))/(rowSums(ems))^2 + 1)
  muY     = log(rowSums(ems)) - sigma2Y/2
  
  logY = log(y)
  zY   = (logY - muY)/sqrt(sigma2Y)
  
  sum(dnorm(zY, log=TRUE))
}

log.normal.mixture.knownRho <- function(delta1, rho1, y, x1, x2){
  
  n1  = ncol(x1)
  n2  = ncol(x2)
  
  if(n1 < 2 || n2 < 2){
    stop("at least two samples are required for x1 and x2\n")
  }
  
  logX1 = log(x1)
  logX2 = log(x2)
  
  mu1 = apply(logX1, 1, mean)
  mu2 = apply(logX2, 1, mean)
  
  
  mu1 = apply(logX1, 1, mean) + log(rho1*delta1)
  mu2 = apply(logX2, 1, mean) + log((1 - rho1)*delta1)
  
  sd1 = apply(logX1, 1, sd)*sqrt((n1 - 1)/n1)
  sd2 = apply(logX2, 1, sd)*sqrt((n2 - 1)/n2)
  
  muX     = cbind(mu1, mu2)
  sigma2X = cbind(sd1*sd1, sd2*sd2)
  
  ems     = exp(muX + sigma2X/2)
  sigma2Y = log(rowSums(ems*ems*(exp(sigma2X) -1))/(rowSums(ems))^2 + 1)
  muY     = log(rowSums(ems)) - sigma2Y/2
  
  logY = log(y)
  zY   = (logY - muY)/sqrt(sigma2Y)

  sum(dnorm(zY, log=TRUE))
}

# ------------------------------------------------------------
# logNormal likelihood for multiple cell types
# ------------------------------------------------------------

log.normal.mixture <- function(rhos, y, X, cellType, borrowInfo4sd=TRUE){
  
  if(min(as.numeric(table(cellType))) < 2){
    stop("at least two samples are required each cell type\n")
  }

  fun1 <- function(v1, t1){
    tapply(log(v1), t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(log(v1), t1, var)
  }
  
  muX  = t(apply(X, 1, fun1, t1=cellType))
  varX = t(apply(X, 1, fun2, t1=cellType))
  nnct = tapply(cellType, cellType, length)
  ntot = sum(nnct)
  
  uct = sort(unique(cellType))
  
  if(any(colnames(muX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(colnames(varX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(names(nnct) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(borrowInfo4sd){
    X1 = log(X)
    for(i in 1:length(uct)){
      wwi = which(cellType == uct[i])
      X1[,wwi] = log(X)[,wwi] - muX[,i]
    }
    varXPop = apply(X1, 1, var)
    
    for(i in 1:length(uct)){
      wi = nnct[i]/ntot
      varX[,i] = wi*varX[,i] + (1-wi)*varXPop
    }
  }

  muX  = t(t(muX) + log(rhos))

  ems  = exp(muX + varX/2)
  varY = log(rowSums(ems*ems*(exp(varX) -1))/(rowSums(ems))^2 + 1)
  muY  = log(rowSums(ems)) - varY/2
  
  logY = log(y)
  zY   = (logY - muY)/sqrt(varY)

  -sum(dnorm(zY, log=TRUE))
}


# ------------------------------------------------------------
# weighted least squares estimates
# ------------------------------------------------------------

ls.log.normal <- function(y, X, cellType, borrowInfo4sd=TRUE, total=1.0){
  
  if(min(as.numeric(table(cellType))) < 2){
    stop("at least two samples are required each cell type\n")
  }
  
  fun1 <- function(v1, t1){
    tapply(log(v1), t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(log(v1), t1, var)
  }
  
  muX  = t(apply(X, 1, fun1, t1=cellType))
  varX = t(apply(X, 1, fun2, t1=cellType))
  nnct = tapply(cellType, cellType, length)
  ntot = sum(nnct)
  
  uct = sort(unique(cellType))
  
  if(any(colnames(muX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(colnames(varX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(names(nnct) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(borrowInfo4sd){
    X1 = log(X)
    for(i in 1:length(uct)){
      wwi = which(cellType == uct[i])
      X1[,wwi] = log(X)[,wwi] - muX[,i]
    }
    varXPop = apply(X1, 1, var)
    
    for(i in 1:length(uct)){
      wi = nnct[i]/ntot
      varX[,i] = wi*varX[,i] + (1-wi)*varXPop
    }
  }

  Z   = exp(muX + varX/2)
  lm1 = lm(y ~ Z)
  
  rhos = coef(lm1)[-1]
  rhos[which(rhos < 0.01)] = 0.01
  rhos = total*rhos/sum(rhos)
  
  list(rhos=rhos)
}

wls.log.normal <- function(y, X, cellType, weight.type=1, borrowInfo4sd=TRUE,
total=1.0)
{
  
  if(min(as.numeric(table(cellType))) < 2){
    stop("at least two samples are required each cell type\n")
  }
  
  fun1 <- function(v1, t1){
    tapply(log(v1), t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(log(v1), t1, var)
  }
  
  muX  = t(apply(X, 1, fun1, t1=cellType))
  varX = t(apply(X, 1, fun2, t1=cellType))
  nnct = tapply(cellType, cellType, length)
  ntot = sum(nnct)
  
  uct = sort(unique(cellType))
  
  if(any(colnames(muX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(colnames(varX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(names(nnct) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(borrowInfo4sd){
    X1 = log(X)
    for(i in 1:length(uct)){
      wwi = which(cellType == uct[i])
      X1[,wwi] = log(X)[,wwi] - muX[,i]
    }
    varXPop = apply(X1, 1, var)
    
    for(i in 1:length(uct)){
      wi = nnct[i]/ntot
      varX[,i] = wi*varX[,i] + (1-wi)*varXPop
    }
  }

  Z   = exp(muX + varX/2)
  lm1 = lm(y ~ Z)

  rhos = coef(lm1)[-1]
  rhos[which(rhos < 0.01)] = 0.01
  rhos = total*rhos/sum(rhos)

  rhos0 = rhos

  for(k in 1:100){
    eta   = drop(Z %*% rhos)

    if(weight.type==1){
      wts  = 1/eta^2
    }else if(weight.type==2){
      Zrho     = t(t(Z) * rhos)
      var.logY = log(rowSums(Zrho*Zrho*(exp(varX) -1))/(eta)^2 + 1)
      mu.logY  = log(eta) - var.logY/2
      varY     = (exp(var.logY) - 1) * exp(2*mu.logY + var.logY)
      wts      = 1/varY
    }else{
      stop("invalid value of weight.type\n")
    }
    
    lm1 = lm(y ~ Z, weights=wts)
    
    rhos = coef(lm1)[-1]
    rhos[which(rhos < 0)] = 0
    rhos = total*rhos/sum(rhos)

    if(max(abs(rhos - rhos0)) < 1e-5){ break }
    
    rhos0 = rhos
  }
  
  
  rhos[which(rhos < 0.01)] = 0.01
  rhos = total*rhos/sum(rhos)
  
  list(rhos=rhos, iter=k, delta=max(abs(rhos - rhos0)))
}

# ------------------------------------------------------------
# glm.log.normal
# ------------------------------------------------------------

glm.log.normal <- function(y, X, cellType, borrowInfo4sd=TRUE, total=1.0){
  
  if(min(as.numeric(table(cellType))) < 2){
    stop("at least two samples are required each cell type\n")
  }
  
  fun1 <- function(v1, t1){
    tapply(log(v1), t1, mean)
  }
  
  fun2 <- function(v1, t1){
    tapply(log(v1), t1, var)
  }
  
  muX  = t(apply(X, 1, fun1, t1=cellType))
  varX = t(apply(X, 1, fun2, t1=cellType))
  nnct = tapply(cellType, cellType, length)
  ntot = sum(nnct)
  
  uct = sort(unique(cellType))
  
  if(any(colnames(muX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(colnames(varX) != uct)){
    stop("cell type mismatch\n")
  }
  
  if(any(names(nnct) != uct)){
    stop("cell type mismatch\n")
  }

  if(borrowInfo4sd){
    X1 = log(X)
    for(i in 1:length(uct)){
      wwi = which(cellType == uct[i])
      X1[,wwi] = log(X)[,wwi] - muX[,i]
    }
    varXPop = apply(X1, 1, var)
    
    for(i in 1:length(uct)){
      wi = nnct[i]/ntot
      varX[,i] = wi*varX[,i] + (1-wi)*varXPop
    }
  }
  
  Z   = exp(muX + varX/2)

  lm1 = lm(y ~ Z)
  
  rhos = coef(lm1)[-1]
  rhos[which(rhos < 0.01)] = 0.01
  rhos = total*rhos/sum(rhos)

  rhos0  = rhos
  eta    = drop(Z %*% rhos)
  theta  = log(eta)
  logY   = log(y)
  sigma2 = 2*sqrt(mean(logY - theta)^2 + 1) - 2
  mu     = theta - sigma2/2

  logLik1  = sum(dnorm(logY, mean=mu, sd=sqrt(sigma2), log=T))
  
  likTrace = logLik1
  rhoTrace = rhos[1]

  for(k in 1:100){
    
    rhos0   = rhos
    logLik0 = logLik1

    ymu    = log(y) - mu
    eta    = drop(Z %*% rhos)
    Delta  = diag(1/eta)
    H      = diag(ymu)
    Ip     = diag(1, nrow=nrow(Z))
    
    dd2 = solve(t(Z) %*% (Delta %*% Delta %*% (Ip + H) %*% Z))
    dd1 = t(Z) %*% Delta %*% ymu

    rhos = rhos0 + drop(dd2 %*% dd1)
    rhos[which(rhos < 0.01)] = 0.01
    rhos = total*rhos/sum(rhos)
    
    if(max(abs(rhos - rhos0)) < 1e-5){ break }
    
    eta    = drop(Z %*% rhos)
    theta  = log(eta)
    sigma2 = 2*sqrt(mean((logY - theta)^2) + 1) - 2
    mu     = theta - sigma2/2
    logLik1 = sum(dnorm(logY, mean=mu, sd=sqrt(sigma2), log=T))
    
    kk = 0
    while(logLik1 < logLik0 ){
      kk = kk + 1
      
      if(kk == 11){
        rhos = rhos0
        break
      }
      
      rhos = rhos0 + drop(dd2 %*% dd1) * (1/2^kk)
      rhos[which(rhos < 0.01)] = 0.01
      rhos = total*rhos/sum(rhos)
      
      if(max(abs(rhos - rhos0)) < 1e-5){ break }
      
      eta    = drop(Z %*% rhos)
      theta  = log(eta)
      sigma2 = 2*sqrt(mean((logY - theta)^2) + 1) - 2
      mu     = theta - sigma2/2
      
      
      logLik1 = sum(dnorm(logY, mean=mu, sd=sqrt(sigma2), log=T))

    }
    
    if(max(abs(rhos - rhos0)) < 1e-5){ break }

    rhoTrace = c(rhoTrace, rhos[1])
    likTrace = c(likTrace, logLik1)
  }
  
  list(rhos=rhos, iter=k, delta=max(abs(rhos - rhos0)))
}

