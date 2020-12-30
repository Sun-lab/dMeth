
library(gridExtra)
library(data.table)
library(ggplot2)
library(cowplot)

#---------------------------------------------------------------------
# a final summary plot
#---------------------------------------------------------------------

cancer_types = c("LUAD", "LUSC", "SKCM")

for(ct1 in cancer_types){
  load(sprintf('_figures_%s/%s_box_plot_cor_RMSE.RData', ct1, ct1))
  
  load(sprintf('_figures_%s/%s_plot_list_B.RData', ct1, ct1))
  p31 = plist[[1]]+theme_cowplot()
  p32 = plist[[2]]+theme_cowplot()
  
  load(sprintf('_figures_%s/%s_plot_list_NK.RData', ct1, ct1))
  p41 = plist[[1]]+theme_cowplot()
  p42 = plist[[2]]+theme_cowplot()
  
  pdf(sprintf('_figures_%s/%s_full.pdf', ct1, ct1), width=8, height = 10)
  grid.arrange(grobs = list(p1_cor, p2_RMSE, p31, p32, p41, p42), ncol = 2)
  dev.off()
  
}

sessionInfo()
gc()

quit(save = 'no')



