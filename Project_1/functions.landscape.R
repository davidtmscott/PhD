#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  March 2021
# Description: functions to generate different landscapes

###############################################################################
#landscapes


simulate.clumpy.landscape<- function(n.patches, extent, clump.length, 
                                     n.foci = NULL, clump.size = NULL){
  # Generates a random clumpy 2D landscape, in a extent*extent square
  # clump.size:  mean number of patches per clump.  Landscape is  uncorrelated if this -> 0
  # clump.length: mean spatial size of clumps
  
  if (is.null(n.foci)){
    # Number of foci for the clumps is a random number, no bigger than n.patches, with mean n.patches/clump.size
    n.foci <- rbinom(1, size = n.patches, prob = 1 / clump.size)
    
    if(n.foci == 0) n.foci <- 1 # cannot be < 1
    
    clump.foci.x <- extent * runif(n.foci)
    clump.foci.y <- extent * runif(n.foci)
  }
  else if (is.matrix(n.foci)){
    clump.foci.x <- n.foci[,1]
    clump.foci.y <- n.foci[,2]
  }
  else {
    clump.foci.x <- extent * runif(n.foci)
    clump.foci.y <- extent * runif(n.foci)
  }

  
  # Use multinomial generator to generate the number of patches per clump, to ensure that
  # the total is n.patches
  clump.n.patches<- c(rmultinom(1, n.patches, prob = rep(1 / n.foci, n.foci))) # c() converts it from matrix to vector
  
  x<- c()
  y<- c()
  for (j in 1:n.foci){
    # Each patch in the focus is positioned at an exponentially distributed distance from the 
    # focus, in a random direction
    
    r<- rexp(clump.n.patches[j], 1 / clump.length) # distance from focus
    theta<- 2 * pi * runif(clump.n.patches[j]) # direction from focus
    xx<- r * cos(theta) + clump.foci.x[j]
    yy<- r * sin(theta) + clump.foci.y[j]
    x<- c(x, xx)
    y<- c(y, yy)
  }
  
  # We now impose periodic boundary conditions, to ensure all patches are in the
  # extent x extent square
  x<- x %% extent
  y <- y %% extent
  
  return(list(x = x, y = y))
  
}


#n.patches = 500
#extent = 10
#
#par(mfrow=c(4,3),
#    mai = c(0.2,0.5,0.2,0.2), 
#    oma = c(0, 0, 0, 0))
#for (clump.length in c(0.25, 0.5, 1, 2)){
#  for (n.foci in c(5, 10, 20)){
#    set.seed(157)
#    data = simulate.clumpy.landscape(n.patches, extent, clump.length, n.foci)
#    all = c(data$x, data$y)
#    range = c(0, extent)
#    #plot(x, y, xlim=range, ylim=range,xlab = "", ylab = "")
#    plot(data$x, data$y, ylim=range, xlim=range, xlab = "", ylab = "")
#    #ggplot(data.frame(data), aes(x, y)) + geom_point() + ylim(range) + xlim(range) + xlab("")+ ylab("")
#  }
#}
#
#mtext("n.foci = 5", side=3, outer=TRUE, line=-1.5, at =c(0.2,-3))
#mtext("n.foci = 10", side=3, outer=TRUE, line=-1.5, at =c(0.5,-3))
#mtext("n.foci = 20", side=3, outer=TRUE, line=-1.5, at =c(0.8,-3))
#
#mtext("clump.length = 0.25", side=2, outer=TRUE, line=-1.5, at =c(0.9,-3))
#mtext("clump.length = 0.5", side=2, outer=TRUE, line=-1.5, at =c(0.65,-3))
#mtext("clump.length = 1", side=2, outer=TRUE, line=-1.5, at =c(0.35,-3))
#mtext("clump.length = 2", side=2, outer=TRUE, line=-1.5, at =c(0.1,-3))


################################################################################


# simulate.corridorey.landscape.1
simulate.corridorey.landscape.1 <- function(n.patches, extent, n.foci, l.noise){
  # Places patches on the edges of a Delaunay triangulation of a set of
  # n.foci random points, with a bit of noise (determined by l.noise)
  
  require(RTriangle)
  x.foci <- c(extent*runif(n.foci - 2), 5, 5) # changed to n.foci - 2
  y.foci <- c(extent*runif(n.foci - 2), 0, 10)
  pp<- pslg(matrix(c(x.foci, y.foci), ncol=2, byrow=F))
  tt<- triangulate(pp, Y=T)

  # place each patch on a randomly chosen edge of the triangulation
  x<- rep(0, n.patches)
  y<- rep(0, n.patches)
  n.edge<- dim(tt$E)[1]
  for (i in 1:n.patches){
    j<- sample(1:n.edge, size=1)
    rr<- runif(1)
    v1<- tt$E[j, 1]
    v2<- tt$E[j, 2]
    
    x[i]<- rr * x.foci[v1] + (1-rr)*x.foci[v2] 
    y[i]<- rr * y.foci[v1] + (1-rr)*y.foci[v2]
  }
  
  # Add some noise
  if (l.noise > 0){
    r.noise<- rexp(n.patches, 1/l.noise)
    theta.noise<- 2*pi*runif(n.patches)
    x.noise <- r.noise * cos(theta.noise)
    y.noise<- r.noise *sin(theta.noise)
    
    # Set noise to zero if it would take the patch out of the landscape
    x.temp<- x + x.noise
    x.noise [x.temp < 0 | x.temp > extent] <- 0
    x<- x + x.noise
    y.temp<- y + y.noise
    y.noise [y.temp < 0 | y.temp > extent]<- 0
    y <- y + y.noise
    
  }
  return(list(x=x, y=y))
}


#n.patches = 500
#extent = 10
#
#par(mfrow=c(3,3),
#    mai = c(0.2,0.5,0.2), 
#    oma = c(0, 0, 0))
#for (l.noise in c(0.1,0.5, 1)){
#  for (n.foci in c(6,12,24)){
#    set.seed(157)
#    range = c(0,extent)
#    plot(simulate.corridorey.landscape.1(n.patches, extent, n.foci, l.noise),
#         ylim=range, xlim=range, xlab = "", ylab = "")
#  }
#}
#
#mtext("n.foci = 6", side=3, outer=TRUE, line=-1.5, at =c(0.2,-3))
#mtext("n.foci = 12", side=3, outer=TRUE, line=-1.5, at =c(0.5,-3))
#mtext("n.foci = 24", side=3, outer=TRUE, line=-1.5, at =c(0.85,-3))
#
#mtext("l.noise = 0.1", side=2, outer=TRUE, line=-1.5, at =c(0.85,-3))
#mtext("l.noise = 0.5", side=2, outer=TRUE, line=-1.5, at =c(0.5,-3))
#mtext("l.noise = 1", side=2, outer=TRUE, line=-1.5, at =c(0.19,-3))

#dev.off()


################################################################################

# was simulate.corridorey.landscape.2
simulate.corridorey.landscape.2 <- function(n.patches, extent, n.foci, l.smooth, n.trials){
  # Generates a random corridorey 2D landscape, in a extent*extent square
  # Favours patches along the zero contour of a random (symmetric about zero) smooth landscape.
  # Creates a surplus of such patches, then returns the most highly ranked
  # n.foci: number of foci that the patches avoid.  Uncorrelated if this -> infinity
  # l.smooth: spatial scale for the smoothing kernel.  Uncorrelated if this -> infinity or zero
  # n.trials: number of trial patches (generated at random)
  
  
  # could set where you want foci to be
  
  x.foci <- extent*runif(n.foci)
  y.foci <- extent*runif(n.foci)
  
  # half of foci are "positive", half negatve
  sign.foci <- rep(c(1, -1), c(n.foci/2, n.foci - n.foci/2))
  
  x.trial<- extent*runif(n.trials)
  y.trial<- extent*runif(n.trials)
  
  score<- rep(0, n.trials)
  
  for (i in 1:n.trials){
    #i=1
    rr<- sqrt((x.trial[i]-x.foci)^2 + (y.trial[i]-y.foci)^2)
    score[i]<- (sum(sign.foci * exp(-rr / l.smooth) / l.smooth^2))^2
  }
  
  sorted <- sort(score, index.return=T)
  
  returned <- sorted$ix[1:n.patches]
  
  return(list(x=x.trial[returned], y=y.trial[returned]))
}

#n.patches = 500
#extent = 10
#
#for (i in c(1.5, 3, 6)){
#  n.trials = n.patches * i
#  par(mfrow=c(3,3),
#      mai = c(0.2,0.5,0.2,0.2), 
#      oma = c(0, 0, 0, 0))
#  for (l.smooth in c(0.1, 1, 10)){
#    for (n.foci in c(2, 5, 10)){
#      set.seed(17)
#      range = c(0,extent)
#      plot(simulate.corridorey.landscape.2(n.patches, extent, n.foci, l.smooth, n.trials),
#           xlim=range, ylim = range, xlab ="", ylab = "")
#    }
#  }
#  mtext("n.foci = 2", side=3, outer=TRUE, line=-1.5, at =c(0.2,-3))
#  mtext("n.foci = 5", side=3, outer=TRUE, line=-1.5, at =c(0.5,-3))
#  mtext("n.foci = 10", side=3, outer=TRUE, line=-1.5, at =c(0.85,-3))
#  
#  
#  
#  mtext("l.smooth = 0.1", side=2, outer=TRUE, line=-1.5, at =c(0.85,-3))
#  mtext("l.smooth = 1", side=2, outer=TRUE, line=-1.5, at =c(0.5,-3))
#  mtext("l.smooth = 10", side=2, outer=TRUE, line=-1.5, at =c(0.19,-3))
#  
#}
 






