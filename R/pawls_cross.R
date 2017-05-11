# without search the whole grid of lambda
pawls_cross = function(x, y, penalty1 = c("L1"), penalty2 = c("L1"), lambda1,lambda2, beta0, w0, 
                       delta, maxIter, intercept = TRUE, criterion = "BIC", startBeta = NULL, startW = NULL) {
    penalty1 <- match.arg(penalty1)
    penalty2 <- match.arg(penalty2)
    L1 <- length(lambda1)
    L2 <- length(lambda2)
    m <- dim(x)[2]
    n <- length(y)
    pre2 <- pre1 <- 0
    lmaxIter <- 50
    beta0 <- SetBeta0(beta0)
    w0 <- ifelse(w0 == 1, 0.99, w0)
    iter <- 0
    index1 <- ceiling(L1/2)
    index2 <- ceiling(L2/2)
    while ((pre1 != index1||pre2 != index2) && (iter < lmaxIter)) {
        iter = iter + 1
        pre2 = index2
        pre1 = index1
        ## fix lambda1
        resw = pawls_grid(x, y, penalty1 = penalty1, penalty2 = penalty2, lambda1=lambda1[index1],lambda2=lambda2,
                          beta0=beta0, w0=w0, delta=delta, maxIter=maxIter, 
                          intercept = intercept, startBeta = startBeta, startW = startW)
        ## search for best lambda2
        bic2 = BIC_cross(as.vector(resw$wloss), apply(matrix(resw$w, L2, n) != 1, 1, sum), 
                         apply(matrix(resw$beta, L2, m) != 0, 1, sum),n) 
        index2 =bic2$index
        crit2 = bic2$crit
        ## fix lambda2
        res = pawls_grid(x, y, penalty1 = penalty1, penalty2 = penalty2, lambda1=lambda1, lambda2=lambda2[index2],
                         beta0=beta0, w0=w0, delta=delta, maxIter=maxIter, 
                         intercept = intercept, startBeta = startBeta, startW = startW)
        ## search for best lambda1
        bic1 = BIC_cross(as.vector(res$wloss), apply(matrix(res$w, L1, n) != 1, 1, sum),
                         apply(matrix(res$beta,L1, m) != 0, 1, sum),n)
        index1 = bic1$index 
        crit1 = bic1$crit
    }
    
    list(beta = as.vector(res$beta[1, index1, ]), w = as.vector(res$w[1, index1, ]),
         opt.lambda1=lambda1[index1], opt.lambda2 = lambda2[index2], iter = iter, crit1=crit1, crit2=crit2)
    
}
