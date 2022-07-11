power_sbm_fun <- function(seed,ns,epsilons,rho,d,a,b,nsims) {
  
  toreturns <-   foreach(epses = c(1:length(epsilons))) %dopar% {
      eps <- epsilons[epses]
      toreturns_epses <- foreach(vals = c(1:length(ns))) %dopar% {
          n <- ns[vals]
          m <- n
          B <- diag(a-b,d) + matrix(b,d,d)
          B2 <- B + diag(eps,d)
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
            Xtrue <- rho*t(assignmentvector1) %*% nus_true1
            Ytrue <-  rho* t(assignmentvector2) %*% nus_true2
            P1 <- Xtrue %*%Ipq %*% t(Xtrue)
            P2 <- Ytrue %*% Ipq %*% t(Ytrue)
            A <- generateAdjacencyMatrix(P1)
            C <- generateAdjacencyMatrix(P2)
            Xhat <- ase(A,d)
            Yhat <- ase(C,d)
            rm(A,C,P1,P2)
          
            Q <- align_matrices(Xhat,Yhat,lambda=.5,eps=.1,niter=20)
            Xnew <- Xhat %*% Q
            test <- nonpar.test(Xnew,Yhat,nsims=500)
            return(test)
          }
          return(vl)
          
        } 
      names(toreturns_epses) <- ns
      return(toreturns_epses)
  }
  
  names(toreturns) <- epsilons
  return(toreturns)
}