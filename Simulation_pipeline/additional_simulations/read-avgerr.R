# sigma2 is the variance of normal CpGs
# lambda = sigma_aberrant^2 / sigma_consistent^2

load("avgerr.RData")
noiselevel = dimnames(avgerr)[1]
for(i in 1:6){
  cat("MSE when Noise level = ", noiselevel[[1]][i],"------------\n")
  print(t(avgerr[i,,,1]))
}