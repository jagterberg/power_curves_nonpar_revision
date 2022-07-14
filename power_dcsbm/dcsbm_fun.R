if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}
if(!require(doParallel)) {
  install.packages("doParallel")
  library(doParallel)
}
numcores<- detectCores()
cl <- makeCluster(4,type = "FORK")
registerDoParallel(cl,numcores)

power_dcsbm_fun <-function(seed,ns,epsilons,rho,d,a,b,nsims) {
  toreturns <-   foreach(epses = c(1:length(epsilons))) %dopar% {
    eps <- epsilons[epses]
    toreturns_epses <- foreach(vals = c(1:length(ns))) %dopar% {
      #toreturns[[epses]][[vals]] <- list()
      n <- ns[vals]
      m <- n
      B <- diag(a-b,d) + matrix(b,d,d)
      B <- rho*B
      B2 <- B
      nus1 <- eigen(B)
      nus_true1 <- nus1$vectors %*% diag(abs(nus1$values)^(1/2),d,d)
      nus2 <- eigen(B2)
      nus_true2 <- nus2$vectors %*% diag(abs(nus2$values)^(1/2),d,d)
      Ipq <- diag(c(1,rep(-1,(d-1))),d,d)
      pis <- rep(1/d,d)
      sigma <- 1/2
      p <- 1
      toReturn <- rep(0,nsims)
      vl <- foreach(iter = c(1:nsims)) %dopar% {
        print(paste0("iter = ",iter, " of ",max(nsims)
                     ,", n = ",n," of ",max(ns),
                     ", eps = ",eps," of ",max(epsilons)))
        assignmentvector1 <- rmultinom(n,1,pis)
        assignmentvector2 <- rmultinom(m,1,pis)
        Xtrue <-.5*t(assignmentvector1) %*% nus_true1
        dgcorrects <- runif(m,.5-eps,.5+eps)
        Ytrue <-  diag(dgcorrects) %*% t(assignmentvector2) %*% nus_true2
        P1 <- Xtrue %*%Ipq %*% t(Xtrue)
        P2 <- Ytrue %*% Ipq %*% t(Ytrue)
        A <- generateAdjacencyMatrix(P1)
        C <- generateAdjacencyMatrix(P2)
        Xhat <- ase(A,d)
        Yhat <- ase(C,d)
        alphahat <- sum(upper.tri(A))/choose(n,2)
        betahat <- sum(upper.tri(C))/choose(n,2)
        Xhat <- Xhat/ sqrt(alphahat)
        Yhat <- Yhat/sqrt(betahat)
        rm(A,C,P1,P2)
        Q <- align_matrices(Xhat,Yhat,lambda=.5,eps=.1,niter=20)
        Xnew <- Xhat %*% Q
        #run test:
        test <- nonpar.test(Xnew,Yhat,nsims=100)
        return(test)
        #toreturns[[epses]][[vals]][[iter]] <- test

      }
      return(vl)
    }
    names(toreturns_epses) <- ns
    return(toreturns_epses)
  }
  names(toreturns) <- epsilons
  return(toreturns)
}