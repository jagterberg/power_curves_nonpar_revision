source("../misc_functions/misc.R")
source("../misc_functions/embed_and_align.R")
source("../misc_functions/align.R")

set.seed(472022)
ns <- seq(200,800,100)
epsilons <- seq(0 ,.2,.1)
sparsities <- seq(.25,1,.25)
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
    y <- 1
    for (rho in sparsities) {
      n <- ns[vals]
      m <- n
      B <- diag(a-b,d) + matrix(b,d,d)
      B2 <- B + diag(eps,d)
      B <- rho * B
      B2 <- rho * B2
      nus1 <- eigen(B)
      nus_true1 <- nus1$vectors %*% diag(abs(nus1$values)^(1/2),d,d)
      nus2 <- eigen(B2)
      nus_true2 <- nus2$vectors %*% diag(abs(nus2$values)^(1/2),d,d)
      Ipq <- diag(c(1,rep(-1,(d-1))),d,d)
      pis <- rep(1/d,d)
      sigma <- 1/2
      p <- 1
      toreturns[[epses]][[vals]][[y]] <- list()
      toReturn <- rep(0,nsims)
      for (iter in c(1:nsims)) {
        print(paste0("iter = ",iter, " of ",max(nsims)
                     ,", n = ",n," of ",max(ns),
                     ", eps = ",eps," of ",max(epsilons)))
        assignmentvector1 <- rmultinom(n,1,pis)
        assignmentvector2 <- rmultinom(m,1,pis)
        Xtrue <- .5*t(assignmentvector1) %*% nus_true1
        Ytrue <-  .5*t(assignmentvector2) %*% nus_true2
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
        test <- nonpar.test(Xnew,Yhat,nsims=500)
        toreturns[[epses]][[vals]][[y]][[iter]] <- test
      }
      y <- y+1
    } 
    
    names(toreturns[[epses]][[vals]]) <- sparsities
  } 
  names(toreturns[[epses]]) <- ns
}
names(toreturns) <- epsilons
save(toreturns,file = "power_sparse_sbm_6-20.Rdata")


ns <- seq(200,800,200)
epsilons <- seq(0 ,.2,.1)
sparsities <- seq(.8,1,.1)

load("power_sparse_sbm_6-20.Rdata")
#grab first eps
toreturnsall <- toreturns
toreturns <- toreturnsall[[1]]
rho_n_matrix1 <- matrix(0,nrow=length(sparsities),ncol=length(ns)) 
row.names(rho_n_matrix1) <- sparsities
colnames(rho_n_matrix1) <- ns
for(i in c(1:length(toreturns))) {
  #for each n
  for(j in c(1:length(toreturns[[i]]))) {
    for (q in c(1:length(toreturns[[i]][[j]]))) {
      if(toreturns[[i]][[j]][[q]]$`estimated p-value` <= .05) {
        #count hte number of times it is over .05.
        rho_n_matrix1[j,i] <- rho_n_matrix1[j,i] + 1
      }
    }
  }
}

rho_n_matrix1/100



