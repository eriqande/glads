
# The code here runs a simulation to explore how a range of processes can generate genomic patterns.
# Motivated by the ideas in Sonya and Kristen's NERC and John Fell Fund grants.

# This version uses functions written in Rcpp which are a bit faster

# First, clear workspace and graphics
rm(list=ls())
graphics.off()
library(gids)
library(abind)

# Higher level parameters that the user will set to run the simulation
initial.population.size <- 200 # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- 1000 # how many linked genes we will deal with
n.alleles.per.locus     <- 5   # we will assume that each locus has the same number of alleles to start with
bv.for.alleles          <- t(array(1:5,c(n.alleles.per.locus,n.loci)))/4-0.25
V.e                     <- 0.1 # standard deviation in environmental component of the phenotype

# a function to generate the initial population structure
initial.struct          <- function(n.N1,n.l,n.a.l){ # n.N1 - initial N, n.l - N loci, n.a.l - alleles / locus
  rand.ints           <- sample(1:n.a.l,n.N1*n.l*2,TRUE)
  struct              <- array(rand.ints,c(n.N1,n.l,2)) # an array
  return(struct)
}



# estimate fitness
fitness <- function(z,n,b0,b1,b2,b3,epsilon){
  # so this is a bit of a wierd fitness function, but it should see if things work or not.
  ze <- z[,1]
  ze2 <- ze^2
  w <- round(b0 + b1*ze + b2*ze2 + b3*n + rnorm(n,0,epsilon),0)
  w <- ifelse(w<0,0,w)
  return(w)
}

# mating function
mating <- function(fitn, struct){
  females      <- which(fitn[,2] %in% 1)
  males        <- which(fitn[,2] %in% 2)
  female.fitn  <- fitn[females,]
  female.stru  <- struct[females,,]
  male.fitn    <- fitn[males,]
  mum.stru     <- female.stru[c(rep(1:length(females),times=female.fitn[,3])),,]
  dads         <- sample(males,sum(female.fitn[,3]),TRUE,prob=male.fitn[,3])
  dad.stru     <- struct[c(dads),,]
  return(list(mum.stru,dad.stru))
}

dispersal <- function(pop1,pop2,rate1to2,rate2to1){
  pops.comb <- abind(pop1,pop2,along=1)
  n1 <- dim(pop1)[1]
  n2 <- dim(pop2)[1]
  p1 <- rep(1,n1)
  p2 <- rep(2,n2)
  p1 <- ifelse(runif(n1)<rate1to2,p2,p1)
  p2 <- ifelse(runif(n2)<rate2to1,p1,p2)
  ps <- c(p1,p2)
  p1 <- pops.comb[ps==1,,]
  p2 <- pops.comb[ps==2,,]
  return(list(p1,p2))
}

start.1 <- struct.1 <- start.2 <- struct.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
res.1 <- res.2 <- rep(NA,2000)


for (i in 1:2000) {
  zs.1        <- rcpp_g2p_map(struct.1, dim(struct.1), bv.for.alleles, n.loci, V.e)
  zs.2        <- rcpp_g2p_map(struct.2, dim(struct.2), bv.for.alleles, n.loci, V.e)
  res.1[i]    <- mean(zs.1[,1])
  res.2[i]    <- mean(zs.2[,1])
  n1          <- dim(zs.1)[1]
  n2          <- dim(zs.2)[1]
  xx          <- paste(i,' ',res.1[i],' ',res.2[i],' ',n1,' ',n2,'\n')
  cat(xx)
  #  w <- round(2 + 0.2*z-0.005*z2-0.01*n + rnorm(n,0,1),0)
  fits.1      <- cbind(zs.1,fitness(zs.1,dim(zs.1)[1],2,0.2,-0.005,-0.01,1))
  pairs.1     <- mating(fits.1,struct.1)
  mums.1      <- pairs.1[[1]]
  dads.1      <- pairs.1[[2]]
  struct.1    <- array(NA, dim(mums.1))
  struct.1[,,1] <- rcpp_recombo_segregate(mums.1, dim(mums.1), rep(0.005, dim(mums.1)[2] - 1))
  struct.1[,,2] <- rcpp_recombo_segregate(dads.1, dim(dads.1), rep(0.005, dim(dads.1)[2] - 1))

  fits.2      <- cbind(zs.2,fitness(zs.2,dim(zs.2)[1],1,0.1,-0.002,-0.00001,1))
  pairs.2     <- mating(fits.2,struct.2)
  mums.2      <- pairs.2[[1]]
  dads.2      <- pairs.2[[2]]
  struct.2    <- array(NA, dim(mums.2))
  struct.2[,,1] <- rcpp_recombo_segregate(mums.2, dim(mums.2), rep(0.005, dim(mums.2)[2] - 1))
  struct.2[,,2] <- rcpp_recombo_segregate(dads.2, dim(dads.2), rep(0.005, dim(dads.2)[2] - 1))

  out <- dispersal(struct.1,struct.2,0.05,0.01)
  struct.1 <- out[[1]]
  struct.2 <- out[[2]]
}









x <- 1:2000
plot(x,res.1,xlab='Generation',ylab='Mean phenotype',type='l')
lines(x,res.2,col='red')
