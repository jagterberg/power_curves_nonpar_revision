source("../misc_functions/misc.R")
source("../misc_functions/embed_and_align.R")
source("../misc_functions/align.R")

set.seed(482022)
ns <-  seq(200,1000,100)
epsilons <- seq(0 ,.3,.1)
sparsities <- seq(.05,.2,.05)
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
    print(paste("n = ",n," eps = ",eps))
    m <- n
    B <- diag(a-b,d) + matrix(b,d,d)
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
    for (iter in c(1:nsims)) {
      print(paste0("iter = ",iter))
      assignmentvector1 <- rmultinom(n,1,pis)
      assignmentvector2 <- rmultinom(m,1,pis)
      Xtrue <-t(assignmentvector1) %*% nus_true1
      dgcorrects <- runif(.5-eps,.5+eps,m)
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
      Q <- align_embeddings(Xhat,Yhat,p=1,q=d-p,lambda=.0001,eps=.001,niter=200,loss="kernel")
      Xnew <- Xhat %*% Q
      #run test:
      test <- nonpar.test(Xnew,Yhat,nsims=1000)
      toreturns[[epses]][[vals]][[iter]] <- test
      
    }
  } 
}
save(toreturns,file = "power_dcsbm_4-8.Rdata")








