source("../misc_functions/misc.R")
source("../misc_functions/embed_and_align.R")
source("../misc_functions/align.R")

set.seed(472022)
ns <-  seq(200,1000,100)
epsilons <- seq(0 ,.3,.1)
#ns <- c(100,200)
toreturns <- list()
d <- 3
a <- .4
b <- .8
nsims <- 100
for (epses in c(1:length(epsilons))) {
  eps <- epsilons[epses]
  toreturns[[epses]] <- list()
  for (vals in c(1:length(ns))) {
    toreturns[[epses]][[vals]] <- list()
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
    for (iter in c(1:nsims)) {
      print(paste0("iter = ",iter, " of ",max(nsims)
                   ,", n = ",n," of ",max(ns),
                   ", eps = ",eps," of ",max(epses)))
      assignmentvector1 <- rmultinom(n,1,pis)
      assignmentvector2 <- rmultinom(m,1,pis)
      Xtrue <-t(assignmentvector1) %*% nus_true1
      Ytrue <-  t(assignmentvector2) %*% nus_true2
      P1 <- Xtrue %*%Ipq %*% t(Xtrue)
      P2 <- Ytrue %*% Ipq %*% t(Ytrue)
      A <- generateAdjacencyMatrix(P1)
      C <- generateAdjacencyMatrix(P2)
      Xhat <- ase(A,d)
      Yhat <- ase(C,d)
      rm(A,C,P1,P2)
      
      #Qinit <- get_sign(Xhat,Yhat)
      #Qinit2 <- Qinit[c((p+1):d),c((p+1):d)]
      Q <- align_matrices_cheap(Xhat,Yhat,lambda=.0001,eps=.001,niter=200)
      #Q <- OTP(Xhat[,c((p+1):d)],Yhat[,c((p+1):d)],Qinit=Qinit2,lambda=.0001,eps=.0001,niter=300)
      #Xnew <- Xhat%*% bdiag(Qinit[1,1],Q)
      Xnew <- Xhat %*% Q
      #run test:
      test <- nonpar.test(Xnew,Yhat,nsims=1000)
      toreturns[[epses]][[vals]][[iter]] <- test
    }
    
  } 
  names(toreturns[[epses]]) <- ns
}
names(toreturns) <- epsilons



save(toreturns,file = "power_sbm_5-10.Rdata")
# 
#  
# load("power1_4-6.Rdata")
# 
 # power1 <- toreturns[[1]]
 # power1_n1 <- power1[[1]]  
 # 
 
# power1_n1$`estimated p-value`  

    
   
    