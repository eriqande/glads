
# The code here runs a simulation to explore how a range of processes can generate genomic patterns.
# Motivated by the ideas in Sonya and Kristen's NERC and John Fell Fund grants.

# This version uses functions written in Rcpp which are a bit faster

# First, clear workspace and graphics
rm(list=ls())
graphics.off()
library(gids)

# Higher level parameters that the user will set to run the simulation
initial.population.size <- 200 # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- 100 # how many linked genes we will deal with
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
fitness <- function(z,n){
  # so this is a bit of a wierd fitness function, but it should see if things work or not.
  ze <- z[,1]
  sex <- z[,2]
  #      up <- 4-0.15*log(n)
  #      lo <- 1-0.15*log(n)
  #      up <- if (up<0) up <- 1 else up
  #      lo <- if (lo<0) up <- 0.1 else lo
  #      w <- ifelse(ze>-2 & ze<4,up,lo)
  #      w <- rpois(n,w)
  w <- round(3 + 0.5*z-0.001*n + rnorm(length(ze),0,0.5),0)
  z2 <- z^2
  w <- round(2 + 0.2*z-0.005*z2-0.01*n + rnorm(n,0,1),0)
  w <- ifelse(w<0,0,w)
  return(w)
}

# mating function
mating <- function(fitn, struct){
  females      <- which(fits[,2] %in% 1)
  males        <- which(fits[,2] %in% 2)
  female.fitn  <- fitn[females,]
  female.stru  <- struct[females,,]
  male.fitn    <- fitn[males,]
  mum.stru     <- female.stru[c(rep(1:length(females),times=female.fitn[,3])),,]
  dads         <- sample(males,sum(female.fitn[,3]),TRUE,prob=male.fitn[,3])
  dad.stru     <- struct[c(dads),,]
  return(list(mum.stru,dad.stru))
}


start <- struct     <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
res <- rep(NA,2000)


for (i in 1:2000) {
  zs <- rcpp_g2p_map(struct, dim(struct), bv.for.alleles, n.loci, V.e)

  res[i]     <- mean(zs[,1])
  cat(paste(dim(zs)[1],' '))
  fits       <- cbind(zs,fitness(zs,dim(zs)[1]))
  pairs      <- mating(fits,struct)
  mums       <- pairs[[1]]
  dads       <- pairs[[2]]

  struct     <- array(NA, dim(mums))
  struct[,,1]     <- rcpp_recombo_segregate(mums, dim(mums), rep(0.005, dim(mums)[2] - 1))
  struct[,,2]     <- rcpp_recombo_segregate(dads, dim(dads), rep(0.005, dim(dads)[2] - 1))
}



x <- 1:2000
plot(x,res,xlab='Generation',ylab='Mean phenotype',type='l')
