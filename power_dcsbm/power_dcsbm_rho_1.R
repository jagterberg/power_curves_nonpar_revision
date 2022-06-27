source("../misc_functions/misc.R")
source("../misc_functions/embed_and_align.R")
source("../misc_functions/align.R")
source("dcsbm_fun.R")

set.seed(482022)
ns <- seq(200,800,100) #seq(200)
epsilons <- seq(0 ,.5,.1)
#sparsities <-  #seq(0,.1,.05)
#ns <- c(100,200)
toreturns <- list()
d <- 3
a <- .4
b <- .8
nsims <- 200

toreturns <- power_dcsbm_fun(6272022,ns,epsilons,rho=1,d,a,b,nsims)

# for (epses in c(1:length(epsilons))) {
#   eps <- epsilons[epses]
#   toreturns[[epses]] <- list()
#   for (vals in c(1:length(ns))) {
#     toreturns[[epses]][[vals]] <- list()
#     n <- ns[vals]
#     m <- n
#     B <- diag(a-b,d) + matrix(b,d,d)
#     B2 <- B 
#     nus1 <- eigen(B)
#     nus_true1 <- nus1$vectors %*% diag(abs(nus1$values)^(1/2),d,d)
#     nus2 <- eigen(B2)
#     nus_true2 <- nus2$vectors %*% diag(abs(nus2$values)^(1/2),d,d)
#     Ipq <- diag(c(1,rep(-1,(d-1))),d,d)
#     pis <- rep(1/d,d)
#     sigma <- 1/2
#     p <- 1
#     toReturn <- rep(0,nsims)
#     for (iter in c(1:nsims)) {
#       print(paste0("iter = ",iter, " of ",max(nsims)
#                    ,", n = ",n," of ",max(ns),
#                    ", eps = ",eps," of ",max(epsilons)))
#       assignmentvector1 <- rmultinom(n,1,pis)
#       assignmentvector2 <- rmultinom(m,1,pis)
#       Xtrue <-.5*t(assignmentvector1) %*% nus_true1
#       dgcorrects <- runif(m,.5-eps,.5+eps)
#       Ytrue <-  diag(dgcorrects) %*% t(assignmentvector2) %*% nus_true2
#       P1 <- Xtrue %*%Ipq %*% t(Xtrue)
#       P2 <- Ytrue %*% Ipq %*% t(Ytrue)
#       A <- generateAdjacencyMatrix(P1)
#       C <- generateAdjacencyMatrix(P2)
#       Xhat <- ase(A,d)
#       Yhat <- ase(C,d)
#       alphahat <- sum(upper.tri(A))/choose(n,2)
#       betahat <- sum(upper.tri(C))/choose(n,2)
#       Xhat <- Xhat/ sqrt(alphahat)
#       Yhat <- Yhat/sqrt(betahat)
#       rm(A,C,P1,P2)
#       Q <- align_matrices(Xhat,Yhat,lambda=.5,eps=.1,niter=20)
#       Xnew <- Xhat %*% Q
#       #run test:
#       test <- nonpar.test(Xnew,Yhat,nsims=500)
#       toreturns[[epses]][[vals]][[iter]] <- test
#       
#     }
#   } 
# }
save(toreturns,file = "power_dcsbm_rho1_6-27.Rdata")


#load("power_dcsbm_6-20.Rdata")
#ns <- seq(200,800,200)
#epsilons <- seq(0 ,.2,.1)
#sparsities <- seq(.05,.2,.05)

# eps_n_matrix <- matrix(0,nrow=length(epsilons),ncol=length(ns)) 
# row.names(eps_n_matrix) <- epsilons
# colnames(eps_n_matrix) <- ns
# for(i in c(1:length(toreturns))) {
#   #for each n
#   for(j in c(1:length(toreturns[[i]]))) {
#     for (q in c(1:length(toreturns[[i]][[j]]))) {
#       if(toreturns[[i]][[j]][[q]]$`estimated p-value` <= .05) {
#         #count hte number of times it is over .05.
#         eps_n_matrix[i,j] <- eps_n_matrix[i,j] + 1
#       }
#     }
#   }
# }
# 
# eps_n_matrix/10
# 
# 
