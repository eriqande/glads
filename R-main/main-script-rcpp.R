
# The code here runs a simulation to explore how a range of processes can generate genomic patterns.
# Motivated by the ideas in Sonya and Kristen's NERC and John Fell Fund grants.

# This version uses functions written in Rcpp which are a bit faster

# First, clear workspace and graphics
rm(list=ls())
graphics.off()
library(glads)
library(abind)
library(parallel)
library(progress)
# library(profvis)
library(ggplot2)
library(dplyr)

# Higher level parameters that the user will set to run the simulation
initial.population.size <- 200 # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- 100 # how many linked genes we will deal with
chromo_mb               <- 2e08
loci.pos                <- sort(floor(runif(n.loci, min = 0, max = chromo_mb)))  # randomly sprinkle the loci along a chromo_mb Mb chromosom
n.alleles.per.locus     <- 5   # we will assume that each locus has the same number of alleles to start with
bv.for.alleles          <- t(array(1:5,c(n.alleles.per.locus,n.loci)))/4-0.25
V.e                     <- 0.01 # standard deviation in environmental component of the phenotype
n.gens                  <- 150 # number of generations.

# a function to generate the initial population structure
initial.struct          <- function(n.N1,n.l,n.a.l){ # n.N1 - initial N, n.l - N loci, n.a.l - alleles / locus
  rand.ints             <- sample(1:n.a.l,n.N1*n.l*2,TRUE)
  struct                <- array(rand.ints,c(n.N1,n.l,2)) # an array
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
  zs             <- rcpp_g2p_map(struct, dim(struct), bvs, nloc, ve)
  fits           <- cbind(zs,fitness(zs,dim(zs)[1],b0,b1,b2,b3,epsil))
  pairs          <- mating(fits,struct)
  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))


#  struct.rt[,,1] <- rcpp_recombo_segregate(mums, dim(mums), rep(0.001, dim(mums)[2] - 1))
#  struct.rt[,,2] <- rcpp_recombo_segregate(dads, dim(dads), rep(0.001, dim(dads)[2] - 1))

  struct.rt[,,1] <- rcpp_recombo_segregate_expo(mums, dim(mums), loci.pos, chromo_mb)
  struct.rt[,,2] <- rcpp_recombo_segregate_expo(dads, dim(dads), loci.pos, chromo_mb)
  return(struct.rt)
}

pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)
start.1 <- struct.1 <- start.2 <- struct.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)

res <- list()

start.1 <- struct.1 <- start.2 <- struct.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)


# testing out the expy_recombo
#pos <- sort(as.integer(floor(runif(n.loci, min = 0, max = 2e08))))
#rcpp_recombo_segregate_expo(struct.1, dim(struct.1), pos, 2e08)



res <- list()
res[["0"]] <- list(struct.1, struct.2)

for (i in 1:n.gens){
  x1 <- list(struct.1,bv.for.alleles,n.loci,V.e,2,0.2,-0.005,-0.01,1)
  x2 <- list(struct.2,bv.for.alleles,n.loci,V.e,1.5,0.1,-0.002,-0.005,1)
  x <- list(x1,x2)
  out <- mclapply(x,doer)
  struct.1 <- out[[1]]
  struct.2 <- out[[2]]
  outd <- fast_dispersal(struct.1,struct.2,0.01,0.01)
  struct.1 <- outd[[1]]
  struct.2 <- outd[[2]]
  if (i %% 10==0) res[[paste(i)]] <- list(struct.1,struct.2)
  pb$tick()
}

pb <- progress_bar$new(format = " Calculating summary statistics [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)
dist <- array(NA,c(n.loci,length(res)))
mean.z <- array(NA,c(length(res),2))
var.z <- array(NA,c(length(res),2))
for (i in 1:length(res)){
  dist[,i] <- apply(res[[i]][[1]],2,mean)-apply(res[[i]][[2]],2,mean)
  z1 <- rcpp_g2p_map(res[[i]][[1]], dim(res[[i]][[1]]), bv.for.alleles, n.loci, V.e)
  z2 <- rcpp_g2p_map(res[[i]][[2]], dim(res[[i]][[2]]), bv.for.alleles, n.loci, V.e)
  mean.z[i,1] <- mean(z1[,1])
  mean.z[i,2] <- mean(z2[,1])
  var.z[i,1]  <- var(z1[,1])
  var.z[i,2]  <- var(z2[,1])
  pb$tick()
}

quartz()
par(mfrow=c(2,2))
x <- (1:length(res))*10

plot(x,mean.z[,1],type='l',ylim=range(mean.z),xlab='Generation',ylab='Mean phenotype')
lines(x,mean.z[,2],col='red')
plot(x,var.z[,1],type='l',ylim=range(var.z),xlab='Generation',ylab='Phenotypic variance')
lines(x,var.z[,2],col='red')

goat <- rep(NA,length(res))
for (i in 1:length(res)){
  goat[i] <- mean(dist[,i])
}

plot(1:n.loci,dist[,length(res)],type='l')
mean(dist[,length(res)])*n.loci/2





# then compute Fst at every locus every 10 generations (that were stored)
# (usig the pegas package....takes a ridiculous amount of time...silly!)
fst_df <- lapply(res, function(x) fst_at_loci_with_pos(x[[1]], x[[2]], loci.pos)) %>%
  bind_rows(.id = "generation") %>%
  mutate(generation = as.numeric(generation))


# plot the Fst values across loci on the last iteration
last_gen <- fst_df %>%
  filter(generation == max(generation))

ggplot(last_gen, aes(x = pos, y = Fst)) +
  geom_point() +
  geom_line()



#plot(1:200,as.vector(fst_df[fst_df$generation=='1500',3]))
