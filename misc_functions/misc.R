Rcpp::cppFunction("
                  NumericMatrix generateAdjacencyMatrix(NumericMatrix pMatrix) {
                  
                  int n = pMatrix.cols();
                  NumericMatrix A(n,n);
                  for(int i = 0; i < n; i ++) {
                  for (int j = i + 1; j < n; j++) {
                  A(i,j) = (int)(rand()%100 < (pMatrix(i,j)* 100));
                  A(j,i) = A(i,j);
                  }
                  }
                  return A;
                  }
                  ")
nonpar.test <- function(Xhat,Yhat,nsims = 100,alpha = .05) {
  dist.mat <- get_dist_matrix(Xhat,Yhat)
  i2 <- setdiff( c(1:(nrow(Xhat) + nrow(Yhat))), c(1:nrow(Xhat)) )
  U <- kernel.stat(Xhat,Yhat,dist=dist.mat,i1 = c(1:nrow(Xhat)),i2 = i2)
  #U <- kernel.stat(Xhat, Ynew)
  testresult <- run_perm_test(U,nsims,Xhat,Yhat,dist.mat = dist.mat,alpha=alpha)
  return(testresult)
}


get_dist_matrix <- function(Z1,Z2,sigma = .5) {
  #new_dat <- rbind(Z1,Z2)
  #D1 <- exp(-(as.matrix(stats::dist(new_dat))^2)/(2*sigma^2))
  m <- nrow(Z2)
  n <- nrow(Z1)
  D1 <- exp(-(as.matrix(stats::dist(Z1))^2)/(2*sigma^2))
  D2 <- exp(-(as.matrix(stats::dist(Z2))^2)/(2*sigma^2))
  D3 <- exp(-rect.dist(Z1,Z2)/(2*sigma^2))
  i1 <- c(1:nrow(Z1))
  i2 <- setdiff(c(1:(m+n)),c(1:nrow(Z1)))
  D <- matrix(0,n+m,n+m)
  D[i1,i1] <- D1
  D[i2,i2] <- D2
  D[i1,i2] <- D3
  D[i2,i1] <- t(D3)
  return(D)
}


rect.dist <- function(X,Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  m <- nrow(Y)
  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))
  D <- tmp2 - 2*tmp1 + tmp3
  #D <- exp(-D/(2*(.5^2)))
  return(D)
}


kernel.stat <- function(X,Y,sigma=0.5,dist = NULL,i1=c(1:nrow(X)),
                        i2=c((nrow(X) + 1):(nrow(X)+nrow(Y)))){
  
  n <- nrow(X)
  m <- nrow(Y)
  
  if (is.null(dist)) {
    tmpXX <- sum(exp(-(as.matrix(stats::dist(X))^2)/(2*sigma^2))) - n
    tmpYY <- sum(exp(-(as.matrix(stats::dist(Y))^2)/(2*sigma^2))) - m
    tmpXY <- sum(exp(-(rect.dist(X,Y))/(2*sigma^2)))
    
    tmp <- tmpXX/(n*(n-1)) + tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)
    
    return((m+n)*tmp)
  } else {
    tmpXX <- sum(dist[i1,i1]) - n
    tmpYY <-  sum(dist[i2,i2]) - m
    tmpXY <- sum(dist[i1,i2])
    tmp<- tmpXX /(n*(n-1)) +tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)
    return((m+n)*tmp)
  }
}

run_perm_test <- function(U,nsims,X,Y,dist.mat = NULL,alpha = .05) {
  toReturn <- rep(-1.0,nsims)
  m <- nrow(X)
  n <- nrow(Y)
  for (i in 1:nsims) {
    #cat(i," out of ",nsims,"\r")
    indices_1 <- sample(c(1:(m + n)),size=n,replace = FALSE)
    indices_2 <- setdiff( c(1:(m+n)), indices_1 )
    
    if(is.null(dist.mat)) {
      dist.mat <- get_dist_matrix(X,Y)
    }
    
    Uhat <- kernel.stat(X=X,Y=Y,i1=indices_1,i2=indices_2,dist=dist.mat)
    #if (!p.val) {
    toReturn[i] <- Uhat
    #} else {
    #  if (Uhat > U) {
    #    toReturn[i] <- 1.0
    #  } else {
    #    toReturn[i] <- 0.0
    #  }
    #}
    
    
  }
  
  crit_val <- sort(toReturn)[floor((1-alpha)*length(toReturn))]
  if (U > crit_val) {
    reject <- "reject the null"
  } else {
    reject <- "do not reject the null"
  }
  
  pval <- sum(toReturn > U)/nsims
  
  final_return <- list(toReturn,reject,crit_val,U,pval)
  names(final_return) <- c("permutation_results","reject","critical_value","test statistic","estimated p-value")
  return(final_return)
  
}


