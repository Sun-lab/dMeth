# sigma2 is the variance of normal CpGs
# lambda = sigma_aberrant^2 / sigma_consistent^2

# install.packages("stargazer")
library(stargazer)

load("avgerr.RData")
noiselevel = dimnames(avgerr)[1]

dim(avgerr)
avgerr[1,,,1]

stargazer(t(avgerr[1,-2,,1]), summary=FALSE)

for(i in 1:6){
  cat("MSE when Noise level = ", noiselevel[[1]][i],"------------\n")
  stargazer(t(avgerr[i,-2,,1]), summary=FALSE)
}


sessionInfo()
gc()

quit(save = 'no')
