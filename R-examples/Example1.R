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
library(glads)


#################################
## 	 Example 1: see Table 1    ##
##	   	recombination map,     ##
## 	 	with fitness function  ##
#################################

###################
### Parameters  ###
###################
initial.population.size=400
n.loci=300					#Number of loci
n.alleles.per.locus=20		#Number of alleles
loci.pos=1:n.loci			#Loci position
add.loci=50					#Number of additive loci
fitness.pos <- 125:(125+(add.loci-1))   #Additive loci position
sex.ratio <- 0.5

bvs = t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus, add.loci))) #Breeding value of additive loci
e.v=0.01 					#Stochastic environmental variant
d.v=1 						#Stochastic demographic variant
n.gens=100					#Number of generations

recom.map= c(rep(0.5, 299))	#Recombination map
Linked<-c(60:69, 150:159, 230:239) #Linked loci
recom.map[Linked]<-0.0001	#Recombination map with linked loci

migration.rate=0			#migration rate


#b0 Maximum generated offspring
#b1 Phenotypic optima
#b2 Variance of the fitness curve
#b3 Density-dependent demographic variant
##Parameters of the fitness function
param.w1 <- list(b0 = 6,b1 = 0.25, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
param.w2 <- list(b0 = 6,b1 = 0.75, b2 = 0.5, b3 = 0.005, d.v = d.v, add.loci = add.loci)

##Parameters of the phenotype function
param.z1 <- list(sex.ratio, fitness.pos, bvs, add.loci, e.v)
param.z2 <- param.z1

###################
### Simulations ###
###################
set.seed(2)
start.1 <- start.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
start <- list(start.1, start.1) #Initial populations

example<-evolve(x = start, time = n.gens, type = "additive", recombination = "map", recom.rate = recom.map, param.z = list(param.z1, param.z2), param.w = list(param.w1, param.w2))

###################
###     Fst     ###
###################
##Commputation of Fst values
#We should first convert the output from class 'struct' to class 'loci'

example.loci <- struct2pegas(example)

#Estimating Fst with library 'pegas'
require(pegas)
FST<-Fst(example.loci)

#Figure
plot(loci.pos, FST[,"Fst"], pch=16, bty="l", xlab="Loci position", ylab="Fst", ylim=c(0,1), type="l", xaxs="i", yaxs="i", las=1, lwd=2  )

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
