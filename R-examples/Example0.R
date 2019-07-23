#rm(list = ls()) #delete later
#devtools::build_manual()
#usethis::use_vignette("glads")
#usethis::use_build_ignore("./src/read-macs-output.cpp")
#usethis::use_build_ignore("./src/read-macs-output.o")

#DEFINIR CLASES
#Verificar funciones con el archivo del apendice (estas son todas viejas)
set.seed(1)
#Initial population structure
initial.population.size=12
n.loci=10
n.alleles.per.locus=2

start <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)

#parameter values
sex.ratio = 0.5 # Unbiased sex ratio.
mean.fitness = 2 #offspring by breeding pair
recom.map=rep(0.5, n.loci-1)
n.gens=10 #number of generations

init.sex <- sample(1:2, initial.population.size, replace=T)

start.1 <- start.2 <- start.3 <- start

x <- list(start.1, start.2, start.3)
s <- list(init.sex, init.sex, init.sex)

system.time(pop<-evolve(x, n.gens, "map", recom.map) )
system.time( pop<-evolve1(start.1, sex.ratio, mean.fitness, recom.map, n.gens, init.sex = NULL) )


pop<-evolve.const(x, n.gens, "average",crossover, loci.pos, chromo_mb)


pop<-evolve(x, n.gens, "dynamic", "map", recom.map, mean.fitness = 3)


dim(x[[1]])[2]

dim(start.1)

dim(pop[[1]])


help("evolve")


Nind=200 	#Number of individuals
Nsnp=1000	#Number of SNPs

#Random genetic identity
n1<-matrix(sample(0:1, Nind*Nsnp,replace=TRUE), Nind, Nsnp)
n2<-matrix(sample(0:1, Nind*Nsnp,replace=TRUE), Nind, Nsnp)

#Diploid structure
L<-list(n1, n2)
struct<-array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)), dimnames=list(rownames(L[[1]]),colnames(L[[1]]), 1:2 ))

#Loci position
colnames(struct)<-sort(sample(1:12000000, 1000, replace=F) )

#########################
###  Parameter values ###
#########################

loci.pos<-as.numeric(colnames(struct) ) #Loci position
chromo_mb<-max(loci.pos) 	 #Chromosome size

sex.ratio = 0.5 			 #Sex ratio
mean.fitness = 3 			 #Offspring per breeding pair
crossover <- 9/100000000.0 	 #Average recombination rate (1.5cM/MB)
mutation.rate=1.1*10^-8 	 #Mutation per site per generation
n.gens=10					 #Number of generations
migration.rate=0.075		 #Migration rate
d.e=0.005					 #Demographic effect to avoid



initial.population.size=100
n.loci=300					#Number of loci
n.alleles.per.locus=20		#Number of alleles
loci.pos=1:n.loci			#Loci position
add.loci=50					#Number of additive loci
add.pos=125					#Starting position of additive loci
sex.ratio = 0.5 			 #Sex ratio

fitness.pos <- 125:174

bv.for.alleles = t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus, add.loci))) #Breeding value of additive loci
e.v=0.01 					#Stochastic environmental variant
V.d=1 						#Stochastic demographic variant
n.gens=100
d.v = 1

b0=6
b1=0.25
b2=0.5
b3 <- d.d <- 0.01

bvs <- bv.for.alleles
recom.map= c(rep(0.5, 299))	#Recombination map
Linked<-c(60:69, 150:159, 230:239) #Linked loci
recom.map[Linked]<-0.0001	#Recombination map with linked loci

set.seed(2)
start.1 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
start.2 <- start.1
start.3 <- start.1

start.1[,,] <-1
start.2[,,] <-2
start.3[,,] <-3

init.sex <- sample(1:2, initial.population.size, replace=T)
sex <- list(init.sex,init.sex,init.sex)
start<-list(start.1,start.1,start.1)

pop<-evolve(start, n.gens, "constant", "map", recom.map, init.sex = sex)
pop<-evolve(start, 100, "constant", "average", 0.001, 1:300, 300, mutation = 0.01, init.sex = NULL)

pop<-evolve(start, n.gens, "dynamic", "average", 0.001, 1:300, 300, param.z = list(sex.ratio = 0.5), param.w = list(mean.fitness = 3, d.d = 0.01))


pop<-evolve(start, n.gens, "custom", "map", recom.map, loci.pos, chromo_mb = 300,
            init.sex = NULL, migration = NULL, mutation.rate = NULL,
            param.z = list(sex.ratio, fitness.pos, bv.for.alleles, 50, e.v),
            param.w = list(b0,b1,b2,b3, d.v, add.loci))

system.time(
pop<-evolve(start, n.gens, "additive", "map", recom.map, loci.pos, chromo_mb = 300,
            init.sex = NULL, migration.rate = 0.001, mutation.rate = NULL,
            param.z = list(sex.ratio, fitness.pos, bvs, add.loci, e.v),
            param.w = list(b0,b1,b2,b3, d.v, add.loci))
)


system.time(
  pop<-evolve(start, 20, "constant", "average", 0.001, 1:300, 300, mutation.rate = NULL, migration.rate = migration)
)

system.time(
  pop<-evolve(start, 100, "constant", "average", 0.001, 1:300, 300, mutation.rate = 0.1, migration.rate = 0.1)
)



migration = matrix(0, 3,3)
diag(migration) = 1

migration[1,2]<-0.05
migration[1,3]<-0.05
migration[3,2]<-0

dim(pop[[3]])

help(evolve)
bvs=bv.for.alleles
do.call("phenotypes", append(list(start.1), param.z))

dim(pop[[1]])

help(evolve)
class(start)

phenotypes(start.1, fitness.pos, bv.for.alleles, 50, e.v)

pop.struct=start.1[,fitness.pos,]
bvs=bv.for.alleles
n.loci.t=50
ve=e.v

phenotype.fx <- function(struct, fitness.pos, bvs, n.loci.t, sex.ratio, ve){
  pop.struct <- struct[ , fitness.pos, ]
  temp <- dim(pop.struct)
  mat <- matrix(1:temp[2],temp[1],temp[2],byrow=TRUE)
  loci.n <- array(mat,c(temp[1],temp[2],2)) # array that gives the locus index at each position in pop.struct
  new.n <- (pop.struct-1)*n.loci.t+loci.n
  bvv <- as.vector(bvs) # turn to bvs
  outp <- array(bvv[new.n],c(temp[1],temp[2],2))
  z <- apply(outp,1,sum)-n.loci.t+rnorm(temp[1],0,ve)
  sex <- rbinom(nrow(struct), 1, sex.ratio)+1
  return(cbind(z = z, sex = sex))
}


fitness.fx <- function(z,sex,n,b0,b1,b2,b3, d.v, n.loci){
  ng <- n.loci
  a <- b0
  b <- b1
  c <- b2
  w <- round( a*exp(-((z-b*ng)^2)/(2*(c*ng)^2) ) - b3*n + rnorm(n,0,d.v) , 0)
  w <- ifelse(w<0,0,w)
  return(w)
}


#######
fitness2 <- function(z,n,b0,b1,b2,b3,epsilon){
  # so this is a bit of a wierd fitness function, but it should see if things work or not.
  ze <- z[,1]
  ze2 <- ze^2
  w <- round(b0 + b1*ze + b2*ze2 + b3*n + rnorm(n,0,epsilon),0)
  w <- ifelse(w<0,0,w)
  return(w)
}





mutate<-function(i, mutation) {
  um <- runif(1)
  um <- ifelse(um > mutation, 0, 1)
  i <- if (um == 0) i <- i else i <- ifelse(i==1,2,1) #Change 0 to 2!!
  return(i)
}


Mutation <- function(genotypes, mutation){
  mutate<- new.genotypes<-apply(genotypes, c(1,2), mutate, mutation = mutation)
  return(new.genotypes)
}

dim(struct.rt[,,1, drop = F])



Mutation2 <- function(genotypes, mutation){
  nl <- dim(genotypes)[2]
  n <- dim(genotypes)[1]
  mut <- array(rbinom(nl*n,1,mutation), dim=c(n, nl,1))
  genotypes[mut==1]<-ifelse(genotypes[mut==1]==1, 0,1)


  return(new.genotypes)
}


genotypes=start.1

n=10
nl=5
array(rbinom(nl*n*2,1,0.1), dim=c(n, nl,2))


struct<-start.1
test<-Mutation(struct, 0.1)

dim(test)

## We first create a random population with 100 individuals and 10 loci
N <- 100  # Population size
nl <- 10  # Number of additive loci
na <- 4  # Number of alleles per locus
G <- initial.struct(N,nl,na)

## Additional parameters needed for the computation of phenotypes
bvs <- t(array( seq(0,1, length = na) ,c(na, nl)))
sex.ratio <- 0.5
e.v=0.01

## Now we compute the additive phenotype value of individuals
phen <- phenotype(G, bvs, nl, sex.ratio, e.v)

## Additional parameteres for the fitness values
b0 <- 6
b1 <- 0.25
b2 <- 0.5
b3 <- 0.01
d.v = 1

## The fitness values of individuals are computed as follows:
fit <- fitness(phen, N, b0, b1, b2, b3, d.v, nl)

dat<-struct2pegas(start, start)

fix(dat)





## We first create a random population with 100 individuals and 10 biallelic loci
N <- 100  # Population size
nl <- 10  # Number of additive loci
na <- 2  # Number of alleles per locus
G <- initial.struct(N,nl,na)

struct2pegas(list(G))




#Example
#We will start with a population of 20 individuals with 10 biallelic SNPs
initial.population.size <- 20
n.loci <- 10
n.alleles.per.locus <- 2

set.seed(1) #setting the seed for reproducible random numbers
start1 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)

##################################
# Type of evolution = "constant" #
##################################

# We set a recombination map and the number of generations to be simulated
recom.map <- rep(0.5, n.loci-1) #all loci are independent
n.gens <- 10
pop <- evolve(list(start1), n.gens, "constant", "map", recom.map)

# A similar simulation but with recombination type "average". We need to specify the position of loci and the size of the chromosome (MB)
loci.pos<- sample(80000: 12000000, 10) #random loci positions
chromo_mb<-12000000 	 #chromosome size
crossover <- 1/100000000.0 #average recombination rate (1 cM/MB)
pop <- evolve(list(start1), n.gens, "constant", "average", crossover, loci.pos, chromo_mb)

# We include mutation rate for the simulations of biallelic loci
pop <- evolve(list(start1), n.gens, "constant", "average", crossover, loci.pos, chromo_mb, mutation.rate = 0.0001)

# A second population is included to incorporate the effect of migration
start2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
pop <- evolve(list(start1, start2), n.gens, "constant", "average", crossover, loci.pos, chromo_mb, migration.rate = 0.01)

# The sex of inidivuals will be set for the starting generations.
sex.start1 <- sample(1:2, initial.population.size, replace=T)
sex.start2 <- sample(1:2, initial.population.size, replace=T)
init.sex <- list(sex.start1, sex.start2)
pop <- evolve(list(start1, start2), n.gens, "constant", "average", crossover, loci.pos, chromo_mb, init.sex = init.sex)


##################################
# Type of evolution = "dynamic"  #
##################################

# We need to set the parameters for the computation of phenotypes and for the fitness function of type 'dynamic'
sex.ratio <- 0.5
mean.fitness <- 3 #mean number of offsprings per breeding pair
d.d <- 0.01 #density-dependent demographic effect
set.seed(1) #setting the seed for reproducible random numbers

pop <- evolve(list(start1), n.gens, "dynamic", "map", recom.map, param.z = list(sex.ratio), param.w = list(mean.fitness, d.d))

##################################
# Type of evolution = "additive" #
##################################

# We need to set the parameters for the computation of phenotypes and for the fitness function of type 'additive'
sex.ratio <- 0.5
fitness.pos <- 1:n.loci #all simulated loci are additive
add.loci <- n.loci
bvs <- t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus, add.loci))) # breeding value of additive loci
e.v <- 0.01 					# stochastic environmental variant
param.z <- list(sex.ratio, fitness.pos, bvs, add.loci, e.v)

b0 <- 6 # maximum number of offspring
b1 <- 1 # phenotypic optima
b2 <- 0.75 # variance of the gaussian curve
b3 <- 0.01 # intensity of the density-dependence
d.v = 1 # stochastic demographic variant
param.w <- list(b0,b1,b2,b3, d.v, add.loci)

set.seed(1) #setting the seed for reproducible random numbers
pop<-evolve(list(start1), n.gens, "additive", "map", recom.map, param.z = param.z, param.w = param.w)

##################################
#  Type of evolution = "custom"  #
##################################

# We need to set a custom 'phenotype' and/or 'fitness' function. In this example the custom functions are very similar to the default ones.

phenotype <- function(struct, ...){
  pop.struct <- struct[ , fitness.pos, ]
  temp <- dim(pop.struct)
  mat <- matrix(1:temp[2],temp[1],temp[2],byrow=TRUE)
  loci.n <- array(mat,c(temp[1],temp[2],2)) # array that gives the locus index at each position in pop.struct
  new.n <- (pop.struct-1)*add.loci+loci.n
  bvv <- as.vector(bvs) # turn to bvs
  outp <- array(bvv[new.n],c(temp[1],temp[2],2))
  z <- apply(outp,1,sum)-add.loci+rnorm(temp[1],0,e.v)
  sex <- rbinom(nrow(struct), 1, sex.ratio)+1
  return(cbind(z = z, sex = sex))
}

fitness <- function(z, sex, n, ...){
  ng <- n.loci
  a <- b0
  b <- b1
  c <- b2
  w <- round( a*exp(-((z-b*ng)^2)/(2*(c*ng)^2) ) - b3*n + rnorm(n,0,d.v) , 0)
  w <- ifelse(w<0,0,w)
  return(w)
}

# We have to assign the new functions into the namespace
assignInNamespace("phenotype",phenotype, ns="glads")
assignInNamespace("fitness",fitness, ns="glads")

set.seed(1) #setting the seed for reproducible random numbers
pop<-evolve(list(start1), n.gens, "custom", "map", recom.map,
            param.z =list(sex.ratio = sex.ratio, fitness.pos = fitness.pos, bvs = bvs, add.loci = add.loci, e.v = e.v),
            param.w = list(b0 = b0, b1 = b1, b2 = b2, b3 = b3, d.v = d.v, n.loci = add.loci))




########################################################################################################
pop<-evolve(list(start1), n.gens, "customZ", "map", recom.map,
            param.z = list(sex.ratio = sex.ratio, fitness.pos = fitness.pos, bvs = bvs, add.loci = add.loci, e.v = e.v),
            param.w = param.w)


pop<-evolve(list(start1), n.gens, "customW", "map", recom.map,
            param.z = param.z,
            param.w = list(b0 = b0, b1 = b1, b2 = b2, b3 = b3, d.v = d.v, n.loci = add.loci))

