#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  October 2020
# Description: mean first passage time functions (mFPT)

################################
# mean first passage time
################################

# has target (sink) no source - any given site can be a source, walker can start anywhere.
#target = sources
#source = targets 


#target = c(11,12)
#source = c(1)
#W[source, target] = 0

## unconditional mFPT, continuous time
## subset landscape (exclude the target site)
#w <- W[-targets, -targets]
## kappa, transition rate to any target states 
#k <- rowSums(W[-targets, targets, drop = FALSE])
#H <- diag(k + rowSums(w))
#T_ <- solve((H - w), rep(1, ncol(w)), tol = exp(-255))


repeat{
  
# continuous time unconditional mean first passage time
mFPT_cont <- function(W, targets){
  # unconditional mFPT, continuous time
  # subset landscape (exclude the target site)
  w <- W[-targets, -targets]
  # kappa, transition rate to any target states 
  k <- rowSums(W[-targets, targets, drop = FALSE])
  H <- diag(k + rowSums(w))
  T_ <- solve((H - w), rep(1, ncol(w)), tol = exp(-255))  # exp(-255)
  return(T_)
}


#log(mFPT_cont(W, targets))[c(1,2,3)]

#mFPT=mFPT_cont(W, targets)[sources]

mFPT_cont_approx<- function(W, targets, diagnostic=T){
  # unconditional mFPT, continuous time
  # subset landscape (exclude the target site)
  # This version calculates an approximation which is appropriate when the landscape
  # consists of a well-connected subnetwork which is poorly connected to the targets.
  # It returns both the mean first passage time, and a diagnostic to say whether 
  # this should be deemed a good approximation
  w <- W[-targets, -targets]
  # kappa, transition rate to any target states 
  k <- rowSums(W[-targets, targets, drop = FALSE])
  
  v<- w-diag(rowSums(w))
  
  ee<- eigen(v)$values[2] #subdominant eigenvalue
  kappa<- sum(k)/length(k)
  
  
  T_ <- 1/kappa
  if(diagnostic){return(list(mFPT=T_, diagnostic=-kappa/ee))}else{return(T_)}
}

#mFPT_cont_approx(W, targets, diagnostic=T)$mFPT
#mFPT_cont_approx(W, targets, diagnostic=T)$diagnostic

# mFPT_cont but with a tryCatch wrapper 

mFPT_cont_TRY <- function(W, target) {
  out <- tryCatch(
    {
      mFPT_cont(W, target)
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





# continuous time unonditional mean squared first passage time
mFPT_cont.2 <- function(W, target, mFPT){
  w <- W[-target, -target]
  # kappa, transition rate to any target states 
  k <- rowSums(W[-target, target, drop = FALSE])
  H <- diag(k + rowSums(w))
  T2 <- 2 * solve((H - w), mFPT, tol = exp(-255))
  return(T2)
}


#-------- for variance:

########################################################
## conditional first passage times

#landscape<- as.data.frame(simulate.clumpy.landscape(1000, 10, 100, 0.5, 5))

# continuous first passage probability
mFPP_cont <- function(W, target, source){
  # remove source and target
  w <- W[-c(target, source), -c(target, source)]
  # internal to target
  k <- rowSums(W[-c(target, source), target, drop = FALSE])
  # internal to source
  n <- rowSums(W[-c(target, source), source, drop = FALSE])
  
  J <- diag(k + n + rowSums(w))
  P <- solve((J - w), k, tol = exp(-255))
  return(P) 
}

#mFPP_cont(W, target, source)

# continuous conditional first moment
mCFPT_cont <- function(W, target, source){
  # conditional mFPT, continuous time
  # first moment
  # subset landscape (exclude the target and source)
  w <- W[-c(target, source), -c(target, source)]
  n <- rowSums(W[-c(target, source), source, drop = FALSE])
  # kappa, transition rate to any target states 
  k <- rowSums(W[-c(target, source), target, drop = FALSE])

  #cbind(cbind(n, w), k)
  
  J <- diag(k + n + rowSums(w))
  P <- mFPP_cont(W, target, source)
  T_1 <- solve((J - w), P, tol = exp(-255))
  return(T_1)
}

#mCFPT_cont(W, target, source)

# continuous conditional second moment
mCFPT_cont.2 <- function(W, target, source){
  # conditional mFPT, continuous time
  # second moment
  w <- W[-c(target, source), -c(target, source)]
  n <- rowSums(W[-c(target, source), source, drop = FALSE])
  k <- rowSums(W[-c(target, source), target, drop = FALSE])
  J <- diag(k + n + rowSums(w))
  T_1 <- mCFPT_cont(W, target, source)
  T_2 <- 2 * solve((J - w), T_1, tol = exp(-255))
  return(T_2)
}

# mCFPT_cont.2(W, target, source)

# conditional variance
vFPT <- function(W, target, source){
  P <- mFPP_cont(W, target, source)
  T1 <- mCFPT_cont(W, target, source) / P # mean FPT
  T2 <- mCFPT_cont.2(W, target, source) / P # mean sq FPT
  
  variance <- T2 - (T1)^2
  return(variance)
}

# unconditional with c() for source
#vFPT(W, target, source)

new_mFPT <- function(W, source, internals, target){
  # probablity of jumping first to internal site i from source k
  Wki <- W[source, internals, drop = FALSE]
  # rate of jumping from source to internals
  Uki <- Wki / rowSums(Wki)
  
  P <- mFPP_cont(W, target, source)
  
  # prob of reaching target from source without returning to source
  #Rk <- colSums(Wki %*% diag(P)) / colSums(Wki)
  Rk <- rowSums(Wki %*% diag(P)) / rowSums(Wki)
  
  T_1 <- mCFPT_cont(W, target, source)
  T_2 <- mCFPT_cont.2(W, target, source)
  
  # test against simulations
  T1 <- rowSums((Uki / Rk) %*% diag(T_1)) + (1 / rowSums(Wki))
  
  T2 <- rowSums((Uki / Rk) %*% diag(T_2)) + 
    rowSums(((2 * Uki) %*% diag(T_1)) / (Rk*rowSums(Wki))) + 
    (2 / rowSums(Wki)^2)
  
  V <- T2 - T1^2
  
  return(list(T1 = T1, T2 = T2, V = V))
}

#

break


################################################################################
# test code from now on
# ignore
################################################################################

break

new_mFPT <- function(W, source, internals, target){
  # probablity of jumping first to internal site i from source k
  Wki <- W[source, internals, drop = FALSE]
  # rate of jumping from source to internals
  if (nrow(Wki) == 1){
    Uki <- Wki / rowSums(Wki)
  } else {
    Uki <- diag(1/rowSums(Wki)) %*% Wki
  }
  
  # implement one note 
  #Wki[2,5] / sum(Wki[2,])
  
  P <- mFPP_cont(W, target, source)
  
  # prob of reaching target from source without returning to source
  #Rk <- colSums(Wki %*% diag(P)) / colSums(Wki)
  Rk <- rowSums(Wki %*% diag(P)) / rowSums(Wki)
  
  T_1 <- mCFPT_cont(W, target, source)
  T_2 <- mCFPT_cont.2(W, target, source)
  
  # test against simulations
  if (nrow(Wki) == 1){
    T1 <- rowSums((Uki / Rk) %*% diag(T_1)) + (1 / rowSums(Wki))
  } else {
    T1 <- rowSums((diag(1/Rk) %*% Uki) %*% diag(T_1)) + (1 / rowSums(Wki))
  }
  
  if (nrow(Wki) == 1){
    T2 <- rowSums((Uki / Rk) %*% diag(T_2)) + 
      rowSums(((2 * Uki) %*% diag(T_1)) / (Rk*rowSums(Wki))) + 
      (2 / rowSums(Wki)^2)
  } else {
    T2 <- rowSums((diag(1/Rk) %*% Uki) %*% diag(T_2)) + 
      rowSums(diag(1/(Rk*rowSums(Wki))) %*% ((2 * Uki) %*% diag(T_1))) + 
      (2 / rowSums(Wki)^2) 
  }
  
  V <- T2 - T1^2
  
  return(list(T1 = T1, T2 = T2, V = V))
}



a = new_mFPT_a(W, source, internals, target)

b = new_mFPT_b(W, source, internals, target)

a$T1
b$T1

# dont want anything else below
stop()


mFPT_disc <- function(W, target){
  # calculates mFPT in discrete time
  # must be given transition probabilities not rates
  
  # subset landscape (exclude the target site)
  w <- W[-target, -target]
  # create identity matrix
  I <- diag(1, nrow(w), ncol(w))
  # subtract transition matrix W from identity matrix I 
  M <- I - w 
  # get inverse of M and multiply by column vector of 1's
  T_ <- solve(M, rep(1, ncol(w)))
  return(T_)
}

#mFPT_disc(W, target)


##################################################################################
# test mFPT_cont with a simulation

source("transition.rate.functions.R")
source("landscape.functions.R")
source("source.target.functions.R")


## set up landscape 

extent = 10
n.patches = 500

n.target = 10
n.source = 10
total.patches = n.patches + n.target + n.source


sources = seq(1, n.source)
targets = seq(total.patches-n.target+1, total.patches)

source = sources
target = targets


if (n.source == 1){
  source_coords <- c(0, 0)
} else{
  source_coords <- data.frame(x=rep(0, n.source),
                              #y=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5))
                              y=seq(from = 0, to = 10, length.out = n.source))
} 


if (n.target == 1){
  target_coords <- c(10,10)
} else{
  target_coords <- data.frame(x=rep(10, n.target),
                              #y=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5))
                              y=seq(from = 0, to = 10, length.out = n.target))
} 

# landscape type
# uniform 
#landscape <- data.frame(x = runif(n.patches, 0, extent), 
#                        y = runif(n.patches, 0, extent))

landscape <- simulate.corridorey.landscape.2(n.patches=n.patches, extent=extent, 
                                             n.foci = 10, l.smooth=0.01, n.trials=n.patches*2)


landscape <- source.target(landscape, 
              source_coords, 
              target_coords)

dispersal = 1

W <- trans_probs(landscape, dispersal, 40, total.patches, extent)

# define list of non target states
internals = as.numeric(colnames(W[,-c(sources, targets)]))


##################
# conditional simulation
##################

# W <- matrix(rnorm(100, 1, 0.45), nrow = 10)
# W <- matrix(runif(100, 0, 1), nrow = 10)
# W <- matrix(rexp(100), nrow = 10)
# W <- matrix(rbinom(100, 10, 0.5), nrow = 10)
# 
# # make it symmetric
# W[lower.tri(W)] = t(W)[lower.tri(W)]
# diag(W) <- 0
# W <- as.data.frame(W)
# colnames(W) <- seq(1:10)

#target = n.patches
#source = 1
#state = 1

next_state_cond <- function(W, state, source, target, step){
  # get list of states
  if (step == 1){
    # cant jump straight to target in first step or to other source
    states = as.numeric(colnames(W[, -c(source, target)]))
  } else {
    # can it go back to itself
    states = as.numeric(colnames(W))#[,-c(as.numeric(state))]))
  }
  # to get probabilities:
  # get trans rates to all other states 
  rates = W[state, states]
  # sum all rates except for i = i
  sumd_rates = sum(rates)
  # randomly select next state based on prob of moving to each state
  next_state = sample(states, size = 1, prob = rates / sumd_rates)
  # generate expo random number for time to transition
  d = rexp(n = 1, rate = sumd_rates)
  
  return(c(next_state, d))
}



#non_targets = 2:9
all_means <- c()
all_means_2 <- c()
all_vars <- c()
all_success <- c()
for (init.state in sources){
  print(init.state)
  all_times = c()
  iterations = 100000 # 1000000
  for (i in 1:iterations){
    time = 0
    state = init.state
    step = 1
    while (state %in% internals | step == 1){
      a <- next_state_cond(W, state, sources, targets, step)
      state = a[1]
      time = time + a[2]
      step = step + 1
    }
    # make it conditional on not going to source
    if(state %in% targets){
      all_times = c(all_times, time)
    }
  }
  # calculate average time across all iterations
  mean_time <- mean(all_times) 
  mean_time_2 <- mean(all_times^2)
  variance_time <- var(all_times)
  success <- length(all_times)
  
  all_means = c(all_means, mean_time)
  all_means_2 = c(all_means_2, mean_time_2)
  all_vars = c(all_vars, variance_time)
  all_success <- c(all_success, success/iterations)
  
  print(sd(all_times))
  print(sd(all_times) / sqrt(length(all_times)))
}


length(all_success)
P <- mFPP_cont(W, targets, sources)

plot(landscape)
T1 <- new_mFPT(W, sources, internals, targets)$T1
T2 <- new_mFPT(W, sources, internals, targets)$T2

all_means
T1
all_means - T1
mean(all_means) - mean(T1)

all_means_2
T2
all_means_2 - T2
mean(all_means_2) - mean(T2)


all_vars
T2 - T1^2
all_vars - (T2 - T1^2)
mean(all_vars) - mean(T2 - T1^2)

################################################################################

##################
# un conditional simulation
##################


# _______
# unconditional mFPT, continuous time
# subset landscape (exclude the target site)
w <- W[-target, -target]
# kappa, transition rate to any target states 
k <- rowSums(W[-target, target, drop = FALSE])
H <- diag(k + rowSums(w))
T_ <- solve((H - w), rep(1, ncol(w)), tol = exp(-255))
# ______




# see if it can return to current state
next_state <- function(W, state, target){
  # get list of states
  states = as.numeric(colnames(W)) #[,-c(as.numeric(state))]))
  # to get probabilities:
  # get trans rates to all other states 
  rates = W[state, states]
  # sum all rates except for i = i
  sumd_rates = sum(rates)
  # randomly select next state based on prob of moving to each state
  next_state = sample(states, size = 1, prob = rates / sumd_rates)
  # generate expo random number for time to transition
  d = rexp(n = 1, rate = sumd_rates)
  
  return(c(next_state, d))
}


#non_targets = 2:9
all_means <- c()
all_means_2 <- c()
all_vars <- c()
for (init.state in sources){
  print(init.state)
  all_times = c()
  iterations = 1000 # 10000000
  for (i in 1:iterations){
    print(i)
    time = 0
    state = init.state
    #while (state != targets){
    while (!(state %in% targets)){
      a <- next_state(W, state, targets)
      state = a[1]
      time = time + a[2]
    }
    all_times = c(all_times, time)
  }
  # calculate average time across all iterations
  mean_time <- mean(all_times) 
  mean_time_2 <- mean(all_times^2)
  variance_time <- var(all_times)
  
  all_means = c(all_means, mean_time)
  all_means_2 = c(all_means_2, mean_time_2)
  all_vars = c(all_vars, variance_time)
}

sd(all_times) / sqrt(length(all_times))


log(mFPT_cont(W, targets)[sources])
log(time)


mFPT <- mFPT_cont(W, targets)
mFPT[sources]
all_means
all_means - mFPT[sources]
mean(all_means) - mean(mFPT[sources])


mFPT2 <- mFPT_cont.2(W, targets, mFPT)
mFPT2[sources]
all_means_2
all_means_2 - mFPT2[sources]
mean(all_means_2) - mean(mFPT2[sources])

mFPT2[sources] - mFPT[sources]^2
(mFPT2[sources] - mFPT[sources]^2) - all_vars
mean((mFPT2[sources] - mFPT[sources]^2)) - mean(all_vars)



################################################################################



}