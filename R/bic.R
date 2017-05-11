## This file is for tunning parameter selection

## BIC

BIC_grid = function(wloss, beta, w) {
  
  l2 <- dim(beta)[1]
  l1 <-  dim(beta)[2]
  index2 <- 1
  index1 <- 1
  BIC.max <- 1e+08
  bicPre = BIC.max
  bicTemp <- matrix(0, l2, l1)
  bicTemp2<- matrix(0, l2, l1)
  wdf <- matrix(0, l2, l1)
  bdf <- matrix(0, l2, l1)
  n <- length(w[1, 1, ])
  pro <- 0.5
  
  for (i in 1 : l2) {
    for (j in 1 : l1) {
      wdf[i, j] <- sum(w[i, j, ] != 1 + 0)
      bdf[i, j] <- sum(beta[i, j, ] != 0 + 0)
      bicTemp[i, j] <- log(wloss[i, j]/(n)) + (bdf[i, j] + wdf[i, j]) * log(n)/(n)
      if(wdf[i,j] >= n * pro || bdf[i,j] + wdf[i,j] >= n){
        bicTemp2[i,j]=-BIC.max
      }else{
        bicTemp2[i,j]=bicTemp[i,j]
        if (bicTemp[i, j] <= bicPre) {
          index2 = i
          index1 = j
          bicPre = bicTemp[i, j]
        }
      }
    }
  }
  bicTemp2 <- ifelse(bicTemp2==-BIC.max,max(bicTemp),bicTemp2)
  res = list(beta = beta[index2, index1, ], w = w[index2,index1,], raw.bic = bicTemp, 
             bic=bicTemp2,index2=index2,index1=index1)
}

BIC_cross = function(loss, dfw, dfb, n ) {
  pro <- 0.5
  df <- dfb + dfw
  fitLoss <- log(loss/n)
  coeff <- log(n) / n
  vl = (fitLoss + coeff * df )
  crit <- vl
  BIC.max=max(vl)
  for(i in 1: length(vl))
  {
    if(dfw[i] >= n * pro || dfb[i] + dfw[i] >= n)
      crit[i] <- BIC.max
  }
  if(min(crit)==BIC.max){
    crit <- vl
  }
  list(index=which.min(crit),crit=crit)
}
