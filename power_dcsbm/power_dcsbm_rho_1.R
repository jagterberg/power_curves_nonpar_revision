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


save(toreturns,file = "power_dcsbm_rho1_6-27.Rdata")

ns <- seq(200,800,100) #seq(200)
epsilons <- seq(0 ,.5,.1)
d <- 3
a <- .4
b <- .8
nsims <- 200

load("power_dcsbm_rho1_6-27.Rdata")

eps_n_matrix <- matrix(0,nrow=length(epsilons),ncol=length(ns)) 
row.names(eps_n_matrix) <- epsilons
colnames(eps_n_matrix) <- ns
for(i in c(1:length(toreturns))) {
  #for each n
  for(j in c(1:length(toreturns[[i]]))) {
    for (q in c(1:length(toreturns[[i]][[j]]))) {
      if(toreturns[[i]][[j]][[q]]$`estimated p-value` <= .05) {
        #count hte number of times it is over .05.
        eps_n_matrix[i,j] <- eps_n_matrix[i,j] + 1
      }
    }
  }
}

eps_n_matrix/nsims
# 
# 
