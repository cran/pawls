
pawls_grid = function(x, y, penalty1 = c("L1"), penalty2 = c("L1"), lambda1, lambda2, beta0, w0, 
                      delta, maxIter, intercept = TRUE, startBeta = NULL, startW = NULL) {
  L1 = length(lambda1)
  L2 = length(lambda2)
  m = dim(x)[2]
  n = dim(x)[1]
  if (is.null(startBeta)) {
    startBeta = rep(0, m)
  }
  if (is.null(startW)) {
    startW = rep(1, n)
  }
  res <- .Call("PAWLS_GRID", x, y, penalty1, penalty2, lambda1, lambda2, beta0, w0, delta, maxIter, 
               ifelse(intercept, 1, 0), startBeta = startBeta, startW = startW)
  
  res = list(beta = array(res[[1]], dim = c(L2, L1, m)), w = array(res[[2]], dim = c(L2, L1, n)), 
             wloss = array(res[[3]], dim = c(L2, L1)), loss = array(res[[4]], dim = c(L2, L1)), 
             iter = array(res[[5]], dim = c(L2, L1)))
  
  res
  
}


