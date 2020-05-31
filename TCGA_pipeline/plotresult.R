library(gridExtra)
library(data.table)
library(ggplot2)
library(cowplot)

setwd('~/GitHub/dMeth/TCGA_pipeline')
setwd('./_figures_SKCM')
load('Cor_SKCM.RData')
corplot <- p1
load('err_SKCM.RData')

tempdata <- melt(as.data.table(err_SKCM))
colnames(tempdata) <- c('Methods','RMSE')
tempdata$cellType = rep(rownames(err_SKCM),5)
p2 <- ggplot(tempdata,aes(x=Methods,y=RMSE)) + geom_boxplot() +
  geom_point(size = 5,aes(colour = cellType)) + theme_cowplot() + ggtitle('RMSE for SKCM')

load('plist_B_SKCM.RData')
p31 <- plist[[1]]+theme_cowplot()
p32 <- plist[[2]]+theme_cowplot()

load('plist_NK_SKCM.RData')
p41 <- plist[[1]]+theme_cowplot()
p42 <- plist[[2]]+theme_cowplot()

pdf('SKCM_full.pdf',width=12,height = 16)
grid.arrange(grobs = list(p1,p2,p31,p32,p41,p42),ncol = 2)
dev.off()
