#' solve the Optimal transport problem
#' @description Optimally transport according to wasserstein cost
#' C_{ij} = d(QX_i, Y_j)^2.  this is the usual optimal transport problem
#' except with an allowance of a parameter Q for an orthogonal matrix.
#' Give two matrices of dimensions n x d and m x y respectively
#' @param X an n x d matrix of points where each row is a point
#' @param Y similar, except possibly a different number of points
#' @param Q An orthogonal matrix
#' @param lambda the parameter to penalize in sinkhorn divergence
#' @param eps the tolerance to end the algorithm
#' @return P, the n x m matrix of assignments
#' @export
optimal_transport <- function(X,Y,lambda = .1,eps = .01) {
  
  
  n <- dim(X)[1]
  m <- dim(Y)[1]
  
  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))
  C <- tmp2 - 2*tmp1 + tmp3
  #C <- as.matrix(rect.dist(X,Y))^2
  
  
  r <- 1/n
  c <- 1/m
  P <- exp(-lambda * C)
  u <- rep(0,n)
  while (max(abs(u - rowSums(P))) > eps) {
    u = rowSums(P)
    P <- r*P / u
    v <- colSums(P)
    P <- c*t( t(P)/v)
  }
  
  return(P)
  
}

#' Procrustes
#' @description Function to find the optimal Q.  
#' @param X the vectors X
#' @param Y the vectors Y
#' @import Matrix
#' @export
#' @return the procrustes solution
procrustes <- function(X,Y) {
  vals <- svd(t(X)%*% Y)
   return(vals$u %*% t(vals$v))
}

opq_project <- function(Q,p=1) {
  d  <- dim(Q)[1]
  Qp <- Q[c(1:p),c(1:p)]
  Qq <- Q[c((p+1):d),c((p+1):d)]
  svd1 <- svd(Qp)
  svd2 <- svd(Qq)
  return(bdiag(svd1$u%*%t(svd1$v),svd2$u%*%t(svd2$v)))
}


#get sign matrix
get_sign <- function(Xhat,Yhat) {
  d <- dim(Xhat)[2]
  ds <- list()
  for ( c in 1:d) {
    ds[[c]] <- c(-1,1)
  }
  signs <- expand.grid(ds)
  Qmin_vanilla <- rep(0,nrow(signs))
  for (allsigns in c(1:nrow(signs))) {
    Qmin_vanilla[allsigns] <- kernel.stat(Xhat,Yhat %*% diag(signs[allsigns,]))
  }
  Qinit <- diag(signs[which.min(Qmin_vanilla),])
  return(Qinit)
}

#Optimal transport procrustes
OTP <- function(Xhat,Yhat,Qinit = NULL,lambda=.1,eps=.01,niter=100) {
  #if Q is null, set to the identity
  if(is.null(Qinit)) {
    d <- dim(Xhat)[2]
    Qinit <- diag(1,d,d)
  }
  
  Xhat <- Xhat%*%Qinit
  Q <- Qinit
  
  for (i in c(1:niter)) {
    #first optimally transport Xhat Q and Y
    Pi <- optimal_transport(Xhat %*%Q ,Yhat,lambda=lambda,eps=eps)
    
    #then procrustes Xhat and Pi Yhat
    Q <- procrustes(Xhat , Pi%*% Yhat)
    
  }
  obj.value <- norm(Xhat %*% Q - Pi %*% Yhat,"F")
  return(list(Q = Q,Pi=Pi,obj.value =obj.value))
}

#function to align over sign matrices but only find the smallest initialization
#assumes p = 1, q = d-1
align_matrices_cheap <- function(Xhat,Yhat,lambda=.1,eps=.01,niter=100) {
  #find Qinit
  d <- dim(Xhat)[2]
  ds <- list()
  for (c in 1:d) {
    ds[[c]] <- c(-1,1)
  }
  signs <- expand.grid(ds)
  obj.values <- rep(0,nrow(signs))
  
  for (allsigns in c(1:nrow(signs))) {
      obj.values[allsigns] <- kernel.stat(Xhat %*% diag(signs[allsigns,]),Yhat)
  }

  Qinit <- diag(signs[which.min(obj.values),])
  Qinit2 <- Qinit[c(2:d),c(2:d)]
  
  #now OTP along negative parts:
  toReturn <- OTP(Xhat[,c(2:d)],Yhat[,c(2:d)],Qinit =Qinit2,lambda=lambda,eps=eps,niter=niter)
  return(bdiag(Qinit[1,1],toReturn$Q))
}



#function to initialize at all sign matrices
align_matrices <- function(Xhat,Yhat,lambda=.1,eps=.01,niter=100,toPrint = FALSE,loss="OTP") {
  d <- dim(Xhat)[2]
  ds <- list()
  for (c in 1:d) {
    ds[[c]] <- c(-1,1)
  }
  signs <- expand.grid(ds)
  obj.values1 <- rep(0,nrow(signs))
  obj.values2 <- rep(0,nrow(signs))
  obj.values3 <- rep(0,nrow(signs))
  Q_news <- list()
  Q_news2 <- list()
  for (allsigns in c(1:nrow(signs))) {
    obj.values1[allsigns] <- kernel.stat(Xhat %*% diag(signs[allsigns,]),Yhat)
    Qinit1 <- diag(signs[allsigns,])
    Qinit2 <- Qinit1[c(2:d),c(2:d)]
    Q_news[[allsigns]] <- OTP(Xhat[,c(2:d)],Yhat[,c(2:d)],Qinit =Qinit2,lambda=lambda,eps=eps,niter=niter)
    obj.values2[allsigns] <- kernel.stat(Xhat %*% bdiag(Qinit1[1,1],Q_news[[allsigns]]$Q),Yhat)
    Q_news2[[allsigns]] <- OTP(Xhat,Yhat,Qinit = Qinit1,lambda=lambda,eps=eps,niter=niter)
    Q_news2[[allsigns]] <- opq_project(Q_news2[[allsigns]]$Q)
    obj.values3[allsigns] <- kernel.stat(Xhat %*% Q_news2[[allsigns]],Yhat)
    
  }
  
  if (min(obj.values1) < min(obj.values2) & min(obj.values1) < min(obj.values3)) {
    return(diag(signs[which.min(obj.values1),])) 
  } else if (min(obj.values2) < min(obj.values3)) {
    minwhich <- which.min(obj.values2)
    q1 <- diag(signs[minwhich,])[1,1]
    Q2 <- Q_news[[minwhich]]$Q
    return(
      bdiag(q1,Q2)
    )
  } else {
    minwhich <- which.min(obj.values3)
    return(Q_news2[[minwhich]])
  }

}
