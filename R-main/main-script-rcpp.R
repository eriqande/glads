
# The code here runs a simulation to explore how a range of processes can generate genomic patterns.
# Motivated by the ideas in Sonya and Kristen's NERC and John Fell Fund grants.

# This version uses functions written in Rcpp which are a bit faster

# First, clear workspace and graphics
rm(list=ls())
graphics.off()
library(gids)
library(abind)
library(parallel)
library(progress)
library(tcltk)
library(profvis)

# Higher level parameters that the user will set to run the simulation
initial.population.size <- 200 # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- 1000 # how many linked genes we will deal with
n.alleles.per.locus     <- 5   # we will assume that each locus has the same number of alleles to start with
bv.for.alleles          <- t(array(1:5,c(n.alleles.per.locus,n.loci)))/4-0.25
V.e                     <- 0.1 # standard deviation in environmental component of the phenotype
n.gens                  <- 100 # number of generations.

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
  #p1 <- rep(1,n1)
  #p2 <- rep(2,n2)
  p1 <- ifelse(runif(n1)<rate1to2,2,1)
  p2 <- ifelse(runif(n2)<rate2to1,1,2)
  ps <- c(p1,p2)
  p1 <- pops.comb[ps==1,,]
  p2 <- pops.comb[ps==2,,]
  return(list(p1,p2))
}


doer <- function(x){
  struct <- x[[1]]
  bvs <- x[[2]]
  nloc <- x[[3]]
  ve <- x[[4]]
  b0 <- x[[5]]
  b1 <- x[[6]]
  b2 <- x[[7]]
  b3 <- x[[8]]
  epsil <- x[[9]]
  zs            <- rcpp_g2p_map(struct, dim(struct), bvs, nloc, ve)
  fits           <- cbind(zs,fitness(zs,dim(zs)[1],b0,b1,b2,b3,epsil))
  pairs          <- mating(fits,struct)
  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))
  struct.rt[,,1] <- rcpp_recombo_segregate(mums, dim(mums), rep(0.005, dim(mums)[2] - 1))
  struct.rt[,,2] <- rcpp_recombo_segregate(dads, dim(dads), rep(0.005, dim(dads)[2] - 1))
  return(struct.rt)
}


#pb <- progress_bar$new(format = " Doing its shit [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)
#pb <- tkProgressBar(title = "progress bar", min = 0, max = n.gens, width = 300)

start.1 <- struct.1 <- start.2 <- struct.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)

res <- list()

n.gens <- 50

profvis({
  for (i in 1:n.gens){
    x1 <- list(struct.1,bv.for.alleles,n.loci,V.e,2,0.2,-0.005,-0.01,1)
    x2 <- list(struct.2,bv.for.alleles,n.loci,V.e,1.5,0.1,-0.002,-0.001,1)
    x <- list(x1,x2)
    out <- mclapply(x,doer)
    struct.1 <- out[[1]]
    struct.2 <- out[[2]]
    set.seed(i)
    outd <- dispersal(struct.1,struct.2,0.1,0.1)

    set.seed(i)
    out2 <- fast_dispersal(struct.1,struct.2,0.1,0.1)
    print(c(i, dim(outd[[1]]), dim(outd[[2]]), dim(out2[[1]]), dim(out2[[2]])))
    print(c(i, all(outd[[1]]==out2[[1]]), all(outd[[2]]==out2[[2]])))
    print(i)
    struct.1 <- outd[[1]]
    struct.2 <- outd[[2]]
    if (i %% 10==0) res[[i/10]] <- list(struct.1,struct.2)
    #pb$tick()
  }
})


