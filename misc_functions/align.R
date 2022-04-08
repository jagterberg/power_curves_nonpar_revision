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

#function to initialize at all sign matrices
align_matrices <- function(Xhat,Yhat,lambda=.1,eps=.01,niter=100,toPrint = FALSE,loss="OTP") {
  ds <- list()
  results <- list()
  p <- dim(Xhat)[2]
  
  if (p > 1) {
    for ( c in 1:p) {
      ds[[c]] <- c(-1,1)
    }
    signs <- expand.grid(ds)
    obj.values <- rep(0,nrow(signs))
    for (allsigns in c(1:nrow(signs))) {
      if(toPrint) {
        print(paste("On sign matrix",allsigns," of ",nrow(signs)))
      }
      results[[allsigns]] <- OTP(Xhat,Yhat,Qinit = diag(signs[allsigns,]),
                                 lambda=lambda,eps=eps,niter=niter)
      if(loss == "OTP") {
        obj.value <- results[[allsigns]]['obj.value']
      } else {
        obj.value <- kernel.stat(Xhat %*% results[[allsigns]]$Q,Yhat)
      }
      
    }
    
    if (loss == "OTP") {
      Qfinal <- results[[which.min(obj.value)]]['Q']
      return(Qfinal)
    } else {
      return(results)
    }
    
  } else {
    signs <- c(1,-1)
    results1 <- OTP(Xhat,Yhat,Qinit=1,lambda=.1,eps=.1,niter=5)
    results2 <- OTP(Xhat,Yhat,Qinit = -1,lambda=.1,eps=.1,niter=5)
    if (results1[['obj.value']] < results2[['obj.value']]) {
      Qfinal <- 1
    } else {
      Qfinal <- -1
    }
    return(Qfinal)
  }
  
 
  
}
