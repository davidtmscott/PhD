#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  October 2020
# Description: condatis functions

################################
# condatis 
################################

# block testing code 


# has target and source sites
# must set source sites
# tweak mean first passage time to incorporate source patches
# target and source column vectors with transition matrix already calculated

#landscape = data.frame(x = runif(10, 0, 10),y = runif(10, 0, 10))
#W <- trans_probs(landscape, 1, 40, 10, 10)
#W <- matrix(runif(100),nrow=10)
#target = 10
#source = 1

#sources = c(1)
#Cin <- colSums(W[sources, -c(targets, sources), drop = FALSE])
## transition rates from from states to target states
#Cout <- rowSums(W[-c(targets, sources), targets, drop = FALSE])
## transition rates within free states
#Cfree <- W[-c(targets, sources), -c(targets, sources)]
#c_sums <- rowSums(cbind(Cin, Cfree, Cout))
#
#M <- diag(c_sums) - Cfree
## voltage (going through each non target or source states)
#v <- solve(M, Cin, tol = exp(-255)) 
## conductance (time to reach target states from sources states)
#K <- v %*% Cout

repeat{

# stephen method
condatis_s <- function(W, targets, sources){
  # transition rate from source to free states
  Cin <- colSums(W[sources, -c(targets, sources), drop = FALSE])
  # transition rates from from states to target states
  Cout <- rowSums(W[-c(targets, sources), targets, drop = FALSE])
  # transition rates within free states
  Cfree <- W[-c(targets, sources), -c(targets, sources)]
  c_sums <- rowSums(cbind(Cfree, Cin, Cout))
  
  M <- diag(c_sums) - Cfree
  # voltage (going through each non target or source states)
  v <- solve(M, Cin, tol = exp(-255)) 
  # conductance (time to reach target states from sources states)
  K <- v %*% Cout
  return(K)
}


# add a tryCatch wrapper to condatis_s
condatis_s_TRY <- function(W, targets, sources) {
  out <- tryCatch(
    {
      condatis_s(W, targets, sources)
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}


#condatis_s(W, target, source)
#condatis_s_TRY(W, targets, sources)

break
################################################################################



# dont want anything else below
stop()

# jenny method
condatis_j <- function(W, target, source){
  # calc dim of inputs to W 
  Cin <- colSums(W[source, -c(target, source), drop = FALSE])
  Cout <- rowSums(W[-c(target, source), target, drop = FALSE])
  Cfree <- W[-c(target, source), -c(target, source)]
  
  M0 <- diag(Cin + Cout + rowSums(Cfree)) - Cfree
  w <- Cin - Cout
  v0 <- solve(M0, w, tol = exp(-255))
  #current in / out 
  I0 <- (v0 + 1) * Cout
  I1 <- (1 - v0) * Cin
  cond <- (sum(I0) + sum(I1)) / 2 # 4
  
  cur <- Cfree * outer(v0, v0, "-")
  flo <- apply(cur, 1, function(x) { sum(abs(x)) })
  flo <- (flo) / 2 + I0 + I1
  flor <- rank(flo, ties.method = "random")
  return(cond)
}

condatis_s(W, targets, sources)
condatis_j(W, targets, sources)

}
