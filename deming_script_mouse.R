#script for running deming regression on globe server
#at first, save logFCmatrixregr, cortestsign500 and deming_minimizer as .Rds objects
#save .Rds objects



#load .Rds objects
logFCmatrixregr = readRDS("logFCmatrixregr_mouse.Rds")
cortestsign500 = readRDS("cortestsign500_mouse.Rds")
totalrownamematrix = readRDS("totalrownamematrix_mouse.Rds")

deming_minimizer = function(logFCmatrixregr){
  fn = function(k_no_first){
    k = c()
    k[1] = 1
    k[2:length(colnames(logFCmatrixregr))] = k_no_first
    res = 0
    for (i in 1:(length(colnames(logFCmatrixregr))-1)){
      namei = colnames(logFCmatrixregr)[i]
      for (j in (i + 1):length(colnames(logFCmatrixregr))){
        namej = colnames(logFCmatrixregr)[j]
        if (cortestsign500[namei, namej] != 1){
          next
        }
        totalrownames = totalrownamematrix[[namei]][[namej]]
        ai = logFCmatrixregr[totalrownames, namei]
        aj = logFCmatrixregr[totalrownames, namej]
        res = res + sum(
          (((aj - (k[j]/k[i])*ai)^2)*((ai - (k[i]/k[j])*aj)^2))/
            (((aj - (k[j]/k[i])*ai)^2)+((ai - (k[i]/k[j])*aj)^2)))/length(totalrownames)
      }
    }
    return(res)
  }
  
  kvec = rnorm(length(colnames(logFCmatrixregr)) - 1, 1, 1)
  ptm <- proc.time()
  optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B")
  proc.time() - ptm
  kres = c(1, optimized$par)
  minimum = optimized$value
  bigres = list(kres, minimum)
  names(bigres) = c("coefs", "minimum")
  return(bigres)
}

#get result
bigres = list()
minimums = c()
for (i in 1:10){
  bigres[[i]] = deming_minimizer(logFCmatrixregr)
  minimums = c(minimums, bigres[[i]]$minimum)
  print(paste0("Opa! ", i, "th minimization done."))
}
kres = bigres[[which.min(minimums)]]$coefs

#dump result as .Rds file
saveRDS(kres, file = "deming_result_mouse.Rds")