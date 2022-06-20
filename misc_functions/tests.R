source("../misc_functions/misc.R")
source("../misc_functions/embed_and_align.R")
source("../misc_functions/align.R")
set.seed(6202022)
ns <- 400
n <- ns[1]
epses <- 0
eps <- epses[1]
d <- 3
a <- .4
b <- .8
nsims <- 20
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
  Q <- align_matrices(Xhat,Yhat,lambda=.5,eps=.1,niter=20)
  #Q <- OTP(Xhat[,c((p+1):d)],Yhat[,c((p+1):d)],Qinit=Qinit2,lambda=.0001,eps=.0001,niter=300)
  #Xnew <- Xhat%*% bdiag(Qinit[1,1],Q)
  Xnew <- Xhat %*% Q
  #run test:
  test <- nonpar.test(Xnew,Yhat,nsims=500)
  toReturn[iter] <- test$`estimated p-value`
}
sum(toReturn <= .1)/nsims - .1
