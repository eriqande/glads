# This version is where Eric is trying it out with biallelic markers and
# trying to seed them with neutral variation from the coalescent


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


# simulate some starting variation from the neutral coalescent (takes about a minute)
# Going from around 1000 SNPs in 400 haploid genomes.
chromo_length <- 1e8
coal_start <- coalescent_starting_variation(n = 800, bp = chromo_length, S = 1000)

# Higher level parameters that the user will set to run the simulation
initial.population.size <- 200 # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- ncol(coal_start$mat) # how many linked genes we will deal with
loci.pos                <- coal_start$pos
n.alleles.per.locus     <- 2   # with snps everything is biallelic
V.e                     <- 0.1 # standard deviation in environmental component of the phenotype
n.gens                  <- 60 # number of generations.


### NOW, I need:
# 1. a function to turn coal_start into a a list of two population structs
# 2. set the breeding values for all but a few loci (or just one locus) to zero.
# 3. fix the fitness function so that each individual has a base relative fitness of
#    one, say, and then we can add to that according to the presence of alleles at selected
#    loci.
# 4. Note that the fitness function has to be different for the two different populations!


# This function takes the output of coalescent_starting_variation and puts it into a list of two structs.
# it also sets the pop attribute to 1 or two on those structs, so that we know what type of fitness
# function to put on it.  It also counts up the number of "2" alleles in each pop at each locus
csv2struct <- function(C, n1) {
  S <- C$mat
  L <- ncol(S)
  n2 <- nrow(S)/2 - n1
  if(n2<=1) {
    stop("n1 is too large")
  }
  s1 <- array(NA, dim = c(n1, L, 2))
  s2 <- array(NA, dim = c(n2, L, 2))

  s1[,,1] <- S[1:n1,]
  s1[,,2] <- S[(n1 + 1):(2 * n1),]
  s2[,,1] <- S[(2 * n1 + 1):(2 * n1 + n2),]
  s2[,,2] <- S[(2 * n1 + n2 + 1):(2 * n1 + n2 + n2),]


  list(
    structs = list(s1, s2),
    freqs1 = colSums(s1[,,1] == 2) + colSums(s1[,,2]==2),
    freqs2 = colSums(s2[,,1] == 2) + colSums(s2[,,2]==2)
  )
}


tmp <- csv2struct(coal_start, initial.population.size)
init.struct <- tmp$structs
init.freqs <- tmp[2:3]

# Now, we set the breeding values.  For now, let's say that there are three SNPs that occur in at least 1 copy,
# but no more than 3 copies in population 2, as close as possible to the center of the chromosom,
# that will confer additive efffects on the phenotype.
candis <- which(init.freqs$freqs2 >= 3 & init.freqs$freqs2 <= 5)  # these are the indexes of the  candidate sites
posis <- loci.pos[candis]  # these are their positions
addis <- candis[order(abs(posis - chromo_length/2))[1]] # indexes of our three sites to use

# now, every allele has a breeding value of 0 except for the "2" allele at those three
# loci, which have a value of 1.  So an individual could ultimately build up to a value of 6
# with three of those loci
bv.for.alleles <- matrix(0, ncol = 2, nrow = n.loci)
bv.for.alleles[addis, 2] <- 1



# estimate fitness.  Eric has fully hacked this so that it depends only on z and b0 in a very simple way
# I could set b0 to around 2 for population 2 and 0 for population 1
fitness <- function(z,n,b0,b1,b2,b3,epsilon){
  ze <- z[,1] * b0
  ze2 <- ze + 1
  w <- rmultinom(n = 1, size = n * 2, prob = ze2/sum(ze2))
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
  zs             <- rcpp_g2p_map(struct, dim(struct), bvs, 0, ve)  # we set the number of loci to zero because we don't need to subtract them off
  zs[,1]         <- abs(zs[,1])  # this is a hack here to accommodate Eric's silly fitness model.
  fits           <- cbind(zs,fitness(zs,initial.population.size,b0,b1,b2,b3,epsil))
  pairs          <- mating(fits,struct)
  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))


#  struct.rt[,,1] <- rcpp_recombo_segregate(mums, dim(mums), rep(0.001, dim(mums)[2] - 1))
#  struct.rt[,,2] <- rcpp_recombo_segregate(dads, dim(dads), rep(0.001, dim(dads)[2] - 1))

  struct.rt[,,1] <- rcpp_recombo_segregate_expo(mums, dim(mums), loci.pos, chromo_length)
  struct.rt[,,2] <- rcpp_recombo_segregate_expo(dads, dim(dads), loci.pos, chromo_length)
  return(struct.rt)
}

pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)
struct.1 <- init.struct[[1]]
struct.2 <- init.struct[[2]]
res <- list()



# testing out the expy_recombo
#pos <- sort(as.integer(floor(runif(n.loci, min = 0, max = 2e08))))
#rcpp_recombo_segregate_expo(struct.1, dim(struct.1), pos, 2e08)



res <- list()
res[["0"]] <- list(struct.1, struct.2)

for (i in 1:n.gens){
  x1 <- list(struct.1,bv.for.alleles,n.loci,V.e,0,0.2,-0.005,-0.01,1)
  x2 <- list(struct.2,bv.for.alleles,n.loci,V.e,10,0.1,-0.002,-0.005,1)
  x <- list(x1,x2)
  out <- mclapply(x,doer)
  struct.1 <- out[[1]]
  struct.2 <- out[[2]]
  outd <- fast_dispersal(struct.1,struct.2,0.01,0.01)
  struct.1 <- outd[[1]]
  struct.2 <- outd[[2]]
  #if (i %% 10==0) res[[paste(i)]] <- list(struct.1,struct.2)
  res[[paste(i)]] <- list(struct.1,struct.2)
  pb$tick()
}


#### HERE IS TIM'S STUFF FOR LOOKING AT TRENDS IN THE PHENOTYPE ####
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





#### COMPUTE FST WITH PLINK ####
fst_df <- lapply(res, function(x) fst_via_plink(x, loci.pos)) %>%
  dplyr::bind_rows(.id = "generation") %>%
  mutate(generation = as.numeric(generation))

fst_df2 <- fst_df %>%
  filter(!is.na(FST)) %>%
  group_by(generation) %>%
  mutate(dens = supsmu(POS, FST)$y)

p <- ggplot(fst_df2, aes(x = POS, y = FST, frame = generation)) +
  geom_point(size = 0.05) +
  geom_line(mapping = aes(x = POS, y = 4*dens), colour = "red")


# you gotta do:  devtools::install_github("dgrtwo/gganimate") to get this one.
library(gganimate)

gganimate(p, "chromo_sim.html")

