#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#
#                      Examples for submitted manuscript:
# The many population genetic and demographic routes to islands of genomic divergence
#
# Example 1. Selection: recombination map and fitness function
# Example 2. Neutral evolution: Mutation rate, average recombination, biallelic SNPs,
#                               without fitness function.
# Questions to:
# Claudio S. Quilodrán. Department of Zoology, University of Oxford.
#         Email: claudio.quilodran@zoo.ox.ac.uk; claudio.quilodran@unige.ch
# Eric C. Anderson. Department of Ecology and Evolutionary Biology, University of California.
#         Email: eric.anderson@noaa.gov
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
rm(list=ls())
graphics.off()
library(gids)
library(progress)


#################################
## 	 Example 1: see Table 1    ##
##	   	recombination map,     ##
## 	 	with fitness function  ##
#################################


##################
### Functions  ###
##################

#########################################################
# a function to generate the initial population structure
initial.struct          <- function(n.N1,n.l,n.a.l){ # n.N1 - initial N, n.l - N loci, n.a.l - alleles / locus
  rand.ints             <- sample(1:n.a.l,n.N1*n.l*2,TRUE)
  struct                <- array(rand.ints,c(n.N1,n.l,2)) # an array
  return(struct)
}

#########################################################
# mating function
resample <- function(x, ...) x[sample.int(length(x), ...)]
mating <- function(fitn, struct){
  females      <- which(fitn[,"sex"] %in% 1)
  males        <- which(fitn[,"sex"] %in% 2)
  try(if(length(females)<1 || length(males)<1 ) stop("Population extinction. There is not enough females or males for mating.", call.=F))
  female.fitn  <- fitn[females,]
  female.stru  <- struct[females,,, drop = F]
  male.fitn    <- fitn[males,]
  mum.stru     <- female.stru[c(rep(1:length(females),times=female.fitn[,"fit"])),,, drop = F]
  dads         <- resample(males,sum(female.fitn[,"fit"]),TRUE,prob=male.fitn[,"fit"])
  dad.stru     <- struct[c(dads),,, drop = F]
  return(list(mum.stru,dad.stru))
}

#########################################################
# dispersal
fast_dispersal <- function(pop1,pop2,rate1to2,rate2to1){
  n1 <- dim(pop1)[1]
  n2 <- dim(pop2)[1]
  L <- dim(pop1)[2]
  g <- 2 #diploidia?
  p1 <- ifelse(runif(n1)<rate1to2,2,1)
  p2 <- ifelse(runif(n2)<rate2to1,1,2)

  ret <- rcpp_dispersal_placement(pop1, pop2, dim(pop1), dim(pop2), p1, p2);
  dim(ret[[1]]) <- c(sum(c(p1, p2)==1), L, g)
  dim(ret[[2]]) <- c(sum(c(p1, p2)==2), L, g)

  ret
}


#########################################################
# fitness function
fitness <- function(z,n,b0,b1,b2,b3,epsilon, n.loci){
  ze <- z[,1]
  ng=n.loci
  a=b0
  b=b1
  c=b2
  w <- round( a*exp(-((ze-b*ng)^2)/(2*(c*ng)^2) ) - b3*n + rnorm(n,0,epsilon) , 0)
  w <- ifelse(w<0,0,w)
  return(w)
}


#########################################################
doer.ex1 <- function(x){
  struct <- x[[1]]
  bvs <- x[[2]]
  nloc <- x[[3]]
  ve <- x[[4]]
  b0 <- x[[5]]
  b1 <- x[[6]]
  b2 <- x[[7]]
  b3 <- x[[8]]
  epsil <- x[[9]]
  theta <- x[[10]]
  add.loci<-x[[11]]
  add.pos<-x[[12]]

  struct.fit<-struct[, add.pos:(add.pos+(add.loci-1)),]


  zs             <- rcpp_g2p_map(struct.fit, dim(struct.fit), bvs, add.loci, ve)
  fits           <- cbind(zs,fit=fitness(zs,dim(zs)[1],b0,b1,b2,b3,epsil, add.loci))
  pairs          <- mating(fits,struct)
  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))

  struct.rt[,,1] <- rcpp_recombo_segregate(mums, dim(mums), theta)
  struct.rt[,,2] <- rcpp_recombo_segregate(dads, dim(dads), theta)

  return(struct.rt)
}

#########################################################
main.function.ex1<- function(Nsize,nloci,nalleles,Ve, Vd, time, fitness.param, migration.rate, recom.map, start.1, start.2, bv.for.alleles, add.loci, add.pos){

  initial.population.size <- Nsize # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
  n.loci                  <- nloci # how many linked genes we will deal with

  loci.pos                <- 1:n.loci

  n.alleles.per.locus     <- nalleles
  bv.for.alleles          <- bv.for.alleles

  add.loci				    <- add.loci
  add.pos					<- add.pos

  V.e                     <- Ve # standard deviation in environmental component of the phenotype
  V.d                     <- Vd # standard demographic variation
  n.gens                  <- time # number of generations.

  disp					<- migration.rate #rate of dispersal
  theta					<- recom.map #recombination rate


  struct.1				<- start.1
  struct.2				<- start.2

  param					<- fitness.param

  param[1,2]

  ############################################
  pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)

  res <- list()
  for (i in 1:n.gens){

    x1 <- list(struct.1,bv.for.alleles,n.loci,V.e,param[1,1],param[1,2],param[1,3],param[1,4],V.d, theta, add.loci, add.pos )
    x2 <- list(struct.2,bv.for.alleles,n.loci,V.e,param[2,1],param[2,2],param[2,3],param[2,4],V.d, theta, add.loci, add.pos )

    x <- list(x1,x2)

    out <- lapply(x, doer.ex1)
    struct.1 <- out[[1]]
    struct.2 <- out[[2]]
    outd <- fast_dispersal(struct.1,struct.2,disp,disp)
    struct.1 <- outd[[1]]
    struct.2 <- outd[[2]]

    if (i %% n.gens==0) res<- list(struct.1,struct.2)
    pb$tick()
  }

  return(res)
}


###################
### Parameters  ###
###################
initial.population.size=400
n.loci=300					#Number of loci
n.alleles.per.locus=20		#Number of alleles
loci.pos=1:n.loci			#Loci position
add.loci=50					#Number of additive loci
add.pos=125					#Starting position of additive loci

bv.for.alleles = t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus, add.loci))) #Breeding value of additive loci
V.e=0.01 					#Stochastic environmental variant
V.d=1 						#Stochastic demographic variant
n.gens=100					#Number of generations

#b0 Maximum generated offspring
#b1 Phenotypic optima
#b2 Variance of the fitness curve
#b3 Density-dependent demographic variant
fitness.param<-cbind(b0=c(6,6), b1=c(0.25,0.75), b2=c(0.5,0.5), b3=c(0.01,0.005)) #Parameter of the fitness function
rownames(fitness.param) <-c("pop1", "pop2")


recom.map= c(rep(0.5, 299))	#Recombination map
Linked<-c(60:69, 150:159, 230:239) #Linked loci
recom.map[Linked]<-0.0001	#Recombination map with linked loci

migration.rate=0			#migration rate


###################
### Simulations ###
###################
set.seed(2)
start.1 <- start.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)

example<-main.function.ex1(initial.population.size, n.loci, n.alleles.per.locus,V.e, V.d, n.gens, fitness.param, migration.rate, recom.map, start.1, start.2, bv.for.alleles, add.loci, add.pos)

#Commputation of Fst values
FST<-fst_at_loci_with_pos(example[[1]] , example[[2]], loci.pos)$Fst

#Figure
plot(loci.pos, FST, pch=16, bty="l", xlab="Loci position", ylab="Fst", ylim=c(0,1), type="l", xaxs="i", yaxs="i", las=1, lwd=2  )

######################################################################################################
######################################################################################################

#################################
## 	 Example 2: Mutation rate, ##
##	   average recombination,  ##
##		   biallelic SNPs      ##
## 	 without fitness function  ##
#################################

##################
### Functions  ###
##################

#########################################################
# Biallelic Mutation
mutate<-function(i){
  um<-runif(1)
  um <- ifelse(um>mutation.rate,0,1)
  i <- if (um==0) i <- i else i <- ifelse(i==1,0,1)
  return(i)
}
Mutation <- function(genotypes, mutation.rate){
  new.genotypes<-apply(genotypes, 2, mutate)
  return(new.genotypes)
}


#########################################################
doer.ex2 <- function(x){

  struct <- x[[1]]
  sex.ratio <- x[[2]]
  mean.fitness <- x[[3]]
  loci.pos <-  x[[4]]
  chromo_mb <- x[[5]]
  d.e <- x[[6]]
  crossover <- x[[7]]
  mutation.rate <- x[[8]]

  z <- as.data.frame(cbind(sex = rbinom(nrow(struct), 1, sex.ratio)+1)) #zs
  rownames(z)=rownames(struct)

  n = nrow(struct)
  z$fit <- round( pmax(rpois(n, lambda=mean.fitness) - d.e*n, 0) )

  pairs <- mating(z,struct)

  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))


  struct.rt[,,1] <- rcpp_recombo_segregate_expo(mums, dim(mums), loci.pos, chromo_mb, crossover)
  struct.rt[,,2] <- rcpp_recombo_segregate_expo(dads, dim(dads), loci.pos, chromo_mb, crossover)

  struct.rt[,,1] <- Mutation(struct.rt[,,1, drop = F], mutation.rate)
  struct.rt[,,2] <- Mutation(struct.rt[,,2, drop = F], mutation.rate)

  return(struct.rt)
}


#########################################################
main.function.ex2<- function(start.1, start.2, sex.ratio, mean.fitness, time, migration.rate, loci.pos, chromo_mb, d.e, crossover, mutation.rate){

  n.gens = time

  struct.1 <- start.1
  struct.2 <- start.2

  sex.ratio <- sex.ratio
  mean.fitness <- mean.fitness

  loci.pos <-  loci.pos
  chromo_mb <- chromo_mb

  disp <- migration.rate
  d.e <- d.e

  crossover <- crossover
  mutation.rate <- mutation.rate

  ############################################
  pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)

  res <- list()
  for (i in 1:n.gens){

    x1 <- list(struct.1, sex.ratio, mean.fitness, loci.pos, chromo_mb, d.e, crossover, mutation.rate)
    x2 <- list(struct.2, sex.ratio, mean.fitness, loci.pos, chromo_mb, d.e, crossover, mutation.rate)

    x <- list(x1,x2)

    out <- lapply(x,doer.ex2)
    struct.1 <- out[[1]]
    struct.2 <- out[[2]]

    outd <- fast_dispersal(struct.1,struct.2,disp,disp)
    struct.1 <- outd[[1]]
    struct.2 <- outd[[2]]
    if (i %% n.gens==0) res<- list(struct.1,struct.2)
    pb$tick()
  }

  return(res)
}


########################
### Generating data  ###
########################
Nind=200 	#Number of individuals
Nsnp=1000	#Number of SNPs

#Random genetic identity
set.seed(1)
n1<-matrix(sample(0:1, Nind*Nsnp,replace=TRUE), Nind, Nsnp)
n2<-matrix(sample(0:1, Nind*Nsnp,replace=TRUE), Nind, Nsnp)

#Diploid structure
L<-list(n1, n2)
struct<-array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)), dimnames=list(rownames(L[[1]]),colnames(L[[1]]), 1:2 ))

#Loci position
colnames(struct)<-sort(sample(1:12000000, 1000, replace=F) )

###################
### Parameters  ###
###################
loci.pos<-as.numeric(colnames(struct) ) #Loci position
chromo_mb<-max(loci.pos) 	 #Chromosome size

sex.ratio = 0.5 			 #Sex ratio
mean.fitness = 3 			 #Offspring per breeding pair
crossover <- 1.0/100000000.0 	 #Average recombination rate (cM/MB)
mutation.rate=1.1*10^-8 	 #Mutation per site per generation
n.gens=100					 #Number of generations
migration.rate=0.075		 #Migration rate
d.e=0.005					 #Demographic effect to avoid exponential growth

start.1=struct #genetic identity of population 1 at time 0
start.2=struct #genetic identity of population 2 at time 0

###################
### Simulations ###
###################
example2<-main.function.ex2(start.1, start.2, sex.ratio, mean.fitness, n.gens, migration.rate, loci.pos, chromo_mb, d.e, crossover, mutation.rate)

#Commputation of Fst values
FST<-fst_at_loci_with_pos(example2[[1]] , example2[[2]], loci.pos)$Fst

#Figure
plot(loci.pos, FST, pch=16, bty="l", xlab="Loci position", ylab="Fst", ylim=c(0,1), type="l", xaxs="i", yaxs="i", las=1, lwd=2  )
