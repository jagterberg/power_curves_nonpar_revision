library(irlba)
ase <- function(A,d ) {
  A_svd <- irlba(A,d)
  Xhat <- A_svd$u %*% diag(A_svd$d)^(1/2)
  return(Xhat)
}

align_embeddings <- function(Xhat,Yhat,p,q,lambda=.0001,eps=.001,niter=200,loss="OTP") {
  if (loss == "OTP") {
    res1 <- align_matrices(as.matrix(Xhat[,c(1:p)]),as.matrix(Yhat[,c(1:p)]))
    res2 <- align_matrices(as.matrix(Xhat[,c((p+1):d)]),as.matrix(Yhat[,c((p+1):d)]),
                           niter=niter,eps=eps,lambda=lambda,toPrint = TRUE,loss=loss)
    Qhat <- bdiag(res1,res2$Q)
  } else {
    res2 <- align_matrices(as.matrix(Xhat[,c((p+1):d)]),as.matrix(Yhat[,c((p+1):d)]),
                           niter=niter,eps=eps,lambda=lambda,toPrint = TRUE,loss=loss)
    obj.value1 <- rep(0,length(res2))
    obj.value2 <- rep(0,length(res2))
    for (signs  in c(1:length(res2))) {
      obj.value1[signs] <- kernel.stat(Xhat %*% bdiag(1,res2[[signs]]$Q),Yhat)
      obj.value2[signs] <- kernel.stat(Xhat %*% bdiag(-1,res2[[signs]]$Q),Yhat)
    }
    
    if (min(obj.value1) < min(obj.value2)) {
      return(bdiag(1,res2[[which.min(obj.value1)]]$Q))
    } else {
      return(bdiag(-1,res2[[which.min(obj.value2)]]$Q))
    }
    
    
    
    
  }
  
  
  
  return(Qhat)
}