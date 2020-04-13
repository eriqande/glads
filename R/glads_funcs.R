#' @title Initial Structure
#'
#' @description This function generates an object of class \code{"struct"} with an initial genetic population structure for simulations.
#' @param N Number of individuals in the initial population
#' @param nl Number of simulated loci
#' @param na Number of alleles at each locus. This parameter set the number of alleles for all loci.
#' @details This function returns a three-dimensional array. Rows represent individuals and columns the different loci. Each element of the array is an integer defining the copy of a given allele at a given locus. The third dimension of the array has two layers representing a pair of homologous chromosomes.
#' @return An object of class \code{"struct"} or an array.
#' @seealso \code{\link{evolve}}
#' @examples
#' ##Initial population size of 10 individuals with 5 biallelic loci
#' initial.population.size=10
#' n.loci=5
#' n.alleles.per.locus=2
#' initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
#' @export
initial.struct <- function(N, nl, na) {
  rand.ints           <- sample(1:na, N * nl * 2, replace = TRUE) # alleles at each locus equifrequent
  struct              <- array(rand.ints, c(N, nl, 2)) # a 3-D array
  class(struct) <- c("struct", "array")
  return(struct)
}


resample <- function(x, ...) x[sample.int(length(x), ...)]
#' Mating function
#'
#' A function to select the mating pairs
#' @param z A data.frame with the sex and fitness value of each individual in the population.
#' @param G An object of class \code{"struct"} with the genetic structure of each individual in the population.
#' @details The function is used internally in the function newborns() in order to select the mating pairs. This function returns a three-dimensional array. Rows represent individuals and columns the different loci. Each element of the array is an integer defining the copy of a given allele at a given locus. The third dimension of the array has two layers representing a pair of homologous chromosomes.
#' @return A list of two objects of class \code{"struct"}. The first element represents the genetic structure of reproductive females and the second element the genetic structure of reproductive males in the population.
#' @seealso \code{\link{newborns}} \code{\link{evolve}}
#' @export
#' @keywords internal
mating <- function(z, G){
  females      <- which(z[,"sex"] %in% 1)
  males        <- which(z[,"sex"] %in% 2)
  try(if(length(females)<1 || length(males)<1 ) stop("Population extinction. There is not enough females or males for mating.", call.=F))
  female.z  <- z[females,]
  female.stru  <- G[females,,, drop = F]
  male.z    <- z[males,]
  mum.stru     <- female.stru[c(rep(1:length(females),times=female.z[,"z"])),,, drop = F]
  dads         <- resample(males,sum(female.z[,"z"]),TRUE,prob=male.z[,"z"])
  dad.stru     <- G[c(dads),,, drop = F]
  return(list(mum.stru,dad.stru))
}


#' Migration function
#'
#' A function for the migration of individuals between populations.
#' @param G1 An object of class \code{"struct"} with the genetic structute of the first population.
#' @param G2 An object of class \code{"struct"} with the genetic structute of the second population.
#' @param rate1to2  A numerical value indicating the rate at which individuals in G1 migrate to G2
#' @param rate2to1 A numerical value indicating the rate at which individuals in G2 migrate to G1
#' @details The function is used internally in the function evolve() in order to allow the migration of individuals between populations. This function returns a list of two-dimensional arrays that represents a pair of homologous chromosomes. Rows represent individuals and columns the different loci. Each element of the array is an integer defining the copy of a given allele at a given locus.
#' @return A list of two objects of class \code{"struct"} representing the new genetic structures of populations after migration.
#' @seealso \code{\link{evolve}}
#' @export
#' @keywords internal
migrate <- function(G1,G2,rate1to2,rate2to1){
  n1 <- dim(G1)[1]
  n2 <- dim(G2)[1]
  L <- dim(G1)[2]
  g <- 2
  p1 <- ifelse(runif(n1)<rate1to2,2,1)
  p2 <- ifelse(runif(n2)<rate2to1,1,2)

  ret <- rcpp_dispersal_placement(G1, G2, dim(G1), dim(G2), p1, p2);
  dim(ret[[1]]) <- c(sum(c(p1, p2)==1), L, g)
  dim(ret[[2]]) <- c(sum(c(p1, p2)==2), L, g)

  ret
}

#' Mutation function
#'
#' This function introduces mutations on biallelic genetic structures.
#' @param G An object of class \code{"struct"} with the genetic structure of individuals before gametogenesis.
#' @param mutation.rate Mutation rate of biallelic loci
#' @details This function is used internally in the function newborns() in order to introduce mutations on biallelic genetic structures before meiosis. The genotypes are introduced as an array of integers defining the copy of a given allele at a given locus. Rows represent individuals and columns the different loci.
#' @return An object of class \code{"struct"} or an array.
#' @seealso \code{\link{evolve}} \code{\link{newborns}}
#' @export
#' @keywords internal
mutation <- function(G, mutation.rate){
  n <- dim(G)[1]
  nl <- dim(G)[2]
  nmut<-rbinom(1, (nl*n*2), mutation.rate)
  pos<-sample(1:(nl*n*2), nmut)
  G[pos] <- G[pos] %% 2+1
  return(G)
}


#' Fitness function
#'
#' This function computes a fitness value (\eqn{\omega}) depending on the phenotype (z) of individuals. The relationship between both variables is assumed to be Gaussian.
#' @param z A data.frame with the phenotype value ('z') and the sex of each individual in a population.
#' @param N Population size of the population.
#' @param b0 A numerical value defining the maximum number of offspring that may be generated by a breeding pair.
#' @param b1 A numerical value defining the phenotypic optima of a given population (see Details).
#' @param b2 A numerical value defining the variance of the Gaussian curve.
#' @param b3 A numerical value defining the intensity of the density-dependence on the fitness of individuals in a population of size 'N'.
#' @param d.v A numerical value defining a stochastic demographic variant in the fitness of individuals.
#' @param add.loci An integer with the total number of additive loci participating in the computation of phenotypes.
#' @details This function is used internally in the function evolve() of type 'selection' in order to compute the fitness value of individuals.
#'
#' The value of 'b1' represents the \eqn{z} value expected to produce the maximum number of offspring.
#' The difference in phenotypic optima between the populations drives the strength of 'divergent selection'. Populations exposed to equal phenotypic optima are considered to be under 'concordant selection'.
#'
#' This fitness function (\eqn{\omega}) has the following form:
#' \deqn{\omega = b_0  \exp^{ -\frac{1}{2} \left ( \frac{4z - b_1n_a}{b_2n_a} \right )^2} - b_3N + \varepsilon_d(0, \sigma_d )}
#' Where \eqn{n_a} is equal to 'add.loci', \eqn{N} is the population size and \eqn{\sigma_d} is equal to 'd.v'. The demographic variant \eqn{\varepsilon_d} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'd.v'.
#' @return A vector with the fitness value (\eqn{\omega}) for each individual.
#' @references
#' Quilodrán, C. S., Ruegg, K., Sendell-Price, A. T., Anderson, E., Coulson, T. and Clegg, S.  (2019).
#' The multiple population genetic and demographic routes to islands of genomic divergence.
#' \emph{bioRxiv}.
#' \doi{10.1101/673483}.
#' @seealso \code{\link{evolve}} \code{\link{phenotype}}
#' @examples
#' ## We first create a random population with 100 individuals and 10 loci
#' N <- 100  # Population size
#' nl <- 10  # Number of additive loci
#' na <- 4  # Number of alleles per locus
#' G <- initial.struct(N,nl,na)
#'
#' ## Additional parameters are needed for the computation of phenotypes
#' bvs <- t(array( seq(0,1, length = na) ,c(na, nl)))
#' sex.ratio <- 0.5
#' e.v=0.01
#'
#' ## Now we compute the additive phenotype value of individuals
#' phen <- phenotype(G, bvs, nl, sex.ratio, e.v)
#'
#' ## Additional parameters are needed for the fitness function
#' b0 <- 6
#' b1 <- 0.25
#' b2 <- 0.5
#' b3 <- 0.01
#' d.v = 1
#'
#' ## The fitness values of individuals are computed as follows:
#' fitness(phen, N, b0, b1, b2, b3, d.v, nl)
#' @export
fitness <- function(z,N,b0,b1,b2,b3, d.v, add.loci){
  ze <- z[,"z"]
  ng <- add.loci
  a <- b0
  b <- b1
  c <- b2
  w <- round( a*exp(-((ze-b*ng)^2)/(2*(c*ng)^2) ) - b3*N + rnorm(N,0,d.v) , 0)
  w <- ifelse(w<0,0,w)
  return(w)
}


#' Phenotype function
#'
#' A function for the computation of additive phenotypes
#' @param G An object of class \code{"struct"} with the genetic structure of each individual in the population.
#' @param bvs A matrix with the breeding values of alleles on each loci. The number of rows is equal to the number of additive loci, while the number of columns is equal to the maximum number of alleles in a locus.
#' @param add.loci An integer with the total number of additive loci participating in the computation of phenotypes.
#' @param sex.ratio A numerical value defining the sex ratio of populations
#' @param e.v A numerical value defining a stochastic environmental variant in the computation of phenotypes.
#' @details This function is used internally in the function evolve() of type 'selection' in order to compute the additive phenotype value of individuals.
#'
#' This phenotype function (\eqn{z}) focuses on an additive genetic genotype-phenotype map. The sum of values of alleles at each locus gives a breeding value (\eqn{b_v}) for each individual at a given locus. The sum of breeding values \eqn{bvs} across loci gives a breeding value for the phenotypes (\eqn{z}), which is computed as follows:
#' \deqn{z = \sum_{v =1 }^{n_a} b_v + \varepsilon_e(0, \sigma_e)}
#' Where \eqn{n_a} is equal to 'add.loci' and \eqn{\sigma_e} is equal to 'e.v'. The environmental contribution \eqn{\varepsilon_e} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'e.v'.
#' @return A data.frame with rows equal to the number of individuals and two columns ('sex' and 'z')
#' @references
#' Quilodrán, C. S., Ruegg, K., Sendell-Price, A. T., Anderson, E., Coulson, T. and Clegg, S.  (2019).
#' The multiple population genetic and demographic routes to islands of genomic divergence.
#' \emph{bioRxiv}.
#' \doi{10.1101/673483}.
#' @seealso \code{\link{evolve}} \code{\link{fitness}}
#' @examples
#' ## We first create a random population with 100 individuals and 10 loci
#' N <- 100  # Population size
#' nl <- 10  # Number of additive loci
#' na <- 4  # Number of alleles per locus
#' G <- initial.struct(N,nl,na)
#'
#' ## Additional parameters are needed for the computation of phenotypes
#' bvs <- t(array( seq(0,1, length = na) ,c(na, nl)))
#' sex.ratio <- 0.5
#' e.v=0.01
#'
#' ## Now we compute the additive phenotype value of individuals
#' phenotype(G, bvs, nl, sex.ratio, e.v)
#' @export
phenotype <- function(G, bvs, add.loci, sex.ratio, e.v){
  dims <- dim(G)
  rcpp_g2p_map(G, dims, bvs, add.loci, sex.ratio, e.v)
}



#' Genetic structure of newborns
#'
#' This function generates the genetic structure of the generation of newborns.
#' @param x List of demographic and genetic parameter values (see Details).
#' @inheritParams evolve
#' @details The function is used internally in the function evolve() in order to generate the genetic structure of the newborns of each generation.
#' A three-dimensional array is returned. Rows represent individuals and columns the different loci. Each element of the array is an integer defining the copy of a given allele at a given locus.
#' The third dimension of the array has two layers representing a pair of homologous chromosomes.
#'
#' The following parameter should be included as a list of parameters 'x':
#' \itemize{
#'   \item{struct:} {An object of class \code{"struct"} with the genetic structure of each individual in the population.}
#'   \item{recom.rate:} {A numerical value of the per-base-pair recombination rate (for example, 1e-08) for recombination type 'average' or a vector of nl - 1 elements with the recombination rate between neigbour loci for type 'map' (nl: number of loci).}
#'   \item{init.sex:} {A vector defining the sex of the initial individuals in the populations (1: females and 2: males). The default value is NULL and assign the sex of the first generation randomly.}
#'   \item{mutation.rate:} {Mutation rate of biallelic loci.}
#'   \item{loci.pos:} {A vector with the position of loci. The default value is NULL, but it is required for recombination type 'average'.}
#'   \item{chromo_mb:} {A numerical value with the size of the simulated chromosome (in megabase). The default value is NULL, but it is required for recombination type 'average'.}
#'   \item{param.z:} {A list with the parameter values for the function computing phenotypes.}
#'   \item{param.w:} {A list with the parameter values for the fitness function.}
#'    }
#'
#' Different types of evolution are available for simulations:
#' \itemize{
#'  \item{'constant'} {a constant population size over time. There is no selection, equal sex ratio and each breeding pair generates two offspring. This case represent neutral evolution in which variation are due to recombination, mutation and migration. Variations in the population sizes are also possible due to the effect of migration.}
#'  \item{'dynamic'} {a dynamic population size over time. This type of evolution introduces to type 'constant' an unequal sex ratio, a variable number of offspring and a density-dependent demographic effect to avoid exponential growth. A list with the parameters for the phenotype function (param.z=list(sex.ratio)) and for the fitness function (param.w=list(mean.fitness, d.d)) are required:
#'    \itemize{
#'      \item{sex.ratio:} {a numerical value defining the sex ratio of populations. This value should be included within a list in 'param.z'.}
#'      \item{mean.fitness:} {an integer with the mean number of offsprings generated by breeding pair. This value should be included within a list in 'param.w'.}
#'      \item{d.d:} {numerical value introducing the density-dependent demographic effect to avoid exponential growth. This value should be included within a list in 'param.w'.}
#'
#'    The sex of individuals (1: females and 2: males) are generated by a binomial distribution with probability of being a male equal to the 'sex.ratio' The fitness of individuals (in number of offspring) are obtained from: \eqn{Poisson(\lambda) - N*d.d}, in which \eqn{\lambda} is equal to 'mean.fitness' and \eqn{N} represents the population size.
#'    }
#'   }
#'   \item{'additive'} {additive phenotype evolution. This function introduce quantitative phenotypes and convergent or divergent selection as implemented in Quilodrán et al. (2019). A list with the parameters for the phenotype function (param.z=list(sex.ratio, fitness.pos, bvs, add.loci, e.v)) and for the fitness function (param.w=list(b0,b1,b2,b3, d.v, add.loci)) are required.
#'
#'   The following parameters are needed for the computation of phenotypes (\eqn{z}):
#'    \itemize{
#'      \item{sex.ratio:} {A numerical value defining the sex ratio of populations. This value should be included within a list in 'param.z'.}
#'      \item{fitness.pos:} {A vector with the position of the additive loci participating in the computation of phenotypes. This value should be included within a list in 'param.z'.}
#'      \item{bvs:} {A matrix with the breeding values of alleles on each loci. The number of rows is equal to the number of additive loci, while the number of columns is equal to the maximum number of alleles in a locus. This object should be included within a list in 'param.z'.}
#'      \item{add.loci:} {An integer with the total number of additive loci participating in the computation of phenotypes. This object should be included within a list in 'param.z'.}
#'      \item{e.v:} {A numerical value defining a stochastic environmental variant in the computation of phenotypes. This value should be included within a list in 'param.z'.}
#'
#'      The default phenotype function (\eqn{z}) focus on an additive genetic genotype-phenotype map. Therefore, the sum of values of alleles at each locus gives a breeding value (\eqn{b_v}) for each individual at a given locus. The sum of breeding values \eqn{bvs} across loci gives a breeding value for the phenotypes (\eqn{z}), which is computed as follows:
#'      \deqn{z = \sum_{v =1 }^{n_a} b_v + \varepsilon_e(0, \sigma_e)}
#'      Where \eqn{n_a} is equal to 'add.loci' and \eqn{\sigma_e} is equal to 'e.v'. The environmental contribution \eqn{\varepsilon_e} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'e.v'.
#'      This function returns a data.frame with rows equal to the number of individuals and two columns ('sex' and 'z')
#'    }
#'   The default fitness function (\eqn{\omega}) computes a Gaussian relationship between \eqn{z} and \eqn{\omega}. The following parameters are needed:
#'    \itemize{
#'      \item{b0:} {A numerical value defining the maximum number of offspring generated by breeding pair. This value should be included within a list in 'param.w'.}
#'      \item{b1:} {A numerical value defining the phenotypic optima. In the gaussiam relationship between \eqn{z} and \eqn{\omega}, the value of 'b1' represent the \eqn{z} value expected to produce the maximum number of offspring. This value should be included within a list in 'param.w'. The difference in phenotypic optima between the populations drives the strength of 'divergent selection'. Populations exposed to equal phenotypic optima are considered to be under 'concordant selection'.}
#'      \item{b2:} {A numerical value defining the variance of the Gaussian curve. This value should be included within a list in 'param.z'.}
#'      \item{b3:} {A numerical value defining the intensity of the density-dependence on the fitness of individuals in a population of size \eqn{N}. This value should be included within a list in 'param.w'.}
#'      \item{d.v:} {A numerical value defining a stochastic demographic variant in the fitness of individuals. This value should be included within a list in 'param.w'.}
#'      \item{add.loci:} {An integer with the total number of additive loci participating in the computation of phenotypes. This object should be included within a list in 'param.w'.}
#'
#'      The default fitness function (\eqn{\omega}) has the form:
#'      \deqn{\omega = b_0  \exp^{ -\frac{1}{2} \left ( \frac{4z - b_1n_a}{b_2n_a} \right )^2} - b_3N + \varepsilon_d(0, \sigma_d )}
#'      Where \eqn{n_a} is equal to 'add.loci', \eqn{N} is the population size and \eqn{\sigma_d} is equal to 'd.v'. The demographic variant \eqn{\varepsilon_d} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'd.v'.
#'      This function returns a vector with the fitness value (\eqn{\omega}) of individuals.
#'    }
#'  }
#'  \item{'custom'} {a custom computation of phenotypes and fitness functions. This functions will defined the type of selection fitting particular case studies. Either a function called 'phenotype' or 'fitness' should be introduced into the namespace, with parameter values included in param.z or param.w.
#'  Assuming a per generation time step, the potential number of offspring produced by each individual depends on its phenotype \eqn{\omega = f(z)}, which in turn depends on the individual genotype and on the environment \eqn{z = g(G, E)}.
#'  \eqn{G} is a numeric value determined by an individual’s genotype, representing the genetic value of the genotype. In the case of an additive genetic map, the genetic value of a genotype will be a breeding value. \eqn{E} represents the effect of the environment on phenotypic expression, and this allows to capture the effects of plasticity on phenotypic expression.
#'  The list of parameters ('param.z' and 'param.w') used for the custom 'phenotype' or 'fitness' function should not include variables. In the first case, the variable phenotype values of individuals are already loaded in the working space as 'z', as well as the genetic structure of each individual in a population, included as the object 'struct' (see Examples). }
#' }
#' The recombination between homologous chromosomes are either of type 'map' or 'average'. The first case needs a vector with the recombination rate (\eqn{\rho}) between neigbour loci of length equal to \eqn{nl - 1} (\eqn{nl}: number of loci). The probability of having a crossover (1) or not (0) is uniformly distributed at a rate defined by the value of \eqn{\rho} between loci
#' (i.e. positions with a probability smaller than \eqn{\rho} recombine). The uniform distribution allows each position with the same values of \eqn{\rho} to have an equal chance of crossover across all iterations. There is no recombination between homologous chromosomes when \eqn{\rho = 0}, both loci are completely linked (e.g. within an inversion or situated close to centromeres), while with a value of \eqn{\rho = 0.5}, the recombination rate is completely random (i.e. both loci are very distant on the same chromosome or are located on different chromosomes).
#' A value of \eqn{\rho < 0.5} means the loci are physically linked. The second case, when recombination is of type 'average', a numerical value with the average recombination rate per base pair should be supplied, the loci position and the size of the chromosome in megabase are also required. The crossover points are exponentially distributed as a Poisson process (see Example).
#' @return An object of class \code{"struct"} or an array.
#' @seealso \code{\link{evolve}}
#' @export
#' @keywords internal
newborns <- function(x, recombination, type){
  struct <- x[[1]]
  recom.rate <- x[[2]]
  init.sex <- x[[3]]
  mutation.rate <- x[[4]]
  loci.pos <- x[[5]]
  chromo_mb <- x[[6]]

  param.z <- x[[7]]
  param.w <- x[[8]]
  fun <- x[[9]]

  switch(type,
         constant = {
           mean.fitness <- 2
           n <- nrow(struct)

           if(n %% 2 == 0){
             z <- as.data.frame(cbind(sex =sample(rep(0:1,each=n/2))+1)) #zs
             rownames(z)=rownames(struct)
           } else {
             neven = floor(n) + floor(n) %% 2 -2
             z <- as.data.frame(cbind(sex =c(sample(rep(0:1,each=neven/2)), sample(0:1,1))+1)) #zs
             rownames(z)=rownames(struct)
           }
           if (!is.null(init.sex)){ z <- as.data.frame(cbind(sex = init.sex)); rownames(z)=rownames(struct) }
           z$z <- rep(mean.fitness, n)
         },
         dynamic = {
           sex.ratio <- param.z[[1]]
           mean.fitness <- param.w[[1]]
           d.d <- param.w[[2]]

           n <- nrow(struct)
           z <- as.data.frame(cbind(sex = rbinom(n, 1, sex.ratio)+1)) #zs
           rownames(z)=rownames(struct)

           z$z <- round( pmax(rpois(n, lambda=mean.fitness) - d.d*n, 0) )
         },
         additive = {
           sex.ratio <- param.z[[1]]
           fitness.pos <- param.z[[2]]
           bvs <- param.z[[3]]
           add.loci <- param.z[[4]]
           e.v <- param.z[[5]]

           b0 <- param.w[[1]]
           b1 <- param.w[[2]]
           b2 <- param.w[[3]]
           b3 <- param.w[[4]]
           d.v <- param.w[[5]]
           n.loci <- param.w[[6]]

           struct.fit <- struct[, fitness.pos,]

           zs <- phenotype(struct.fit, bvs, add.loci, sex.ratio, e.v)
           fit <- fitness(zs,dim(zs)[1], b0, b1, b2, b3, d.v, n.loci)
           z <- cbind(sex = zs[,"sex"], z = fit)
          },
         custom = {
           zs <- do.call(fun[1], append(list(struct=struct), param.z))
           fit <- do.call(fun[2], append(list(z = zs[,"z"], sex = zs[,"sex"], n = dim(zs)[1]), param.w))
           z <- cbind(sex = zs[,"sex"], z = fit)
         },
         stop("Invalid type of evolution. Current options are 'constant', 'dynamic', 'additive', and 'custom'")
  )


  pairs <- mating(z,struct)

  mums           <- pairs[[1]]
  dads           <- pairs[[2]]

  if (!is.null(mutation.rate)) {
    mums<- mutation(mums, mutation.rate)
    dads<- mutation(dads, mutation.rate)
  }

  struct.rt      <- array(NA, dim(mums))

  switch(recombination,
         map = {
           struct.rt[,,1] <- rcpp_recombo_segregate(mums, dim(mums), recom.rate)
           struct.rt[,,2] <- rcpp_recombo_segregate(dads, dim(dads), recom.rate)
         },
         average = {
           struct.rt[,,1] <- rcpp_recombo_segregate_expo(mums, dim(mums), loci.pos, chromo_mb, recom.rate)
           struct.rt[,,2] <- rcpp_recombo_segregate_expo(dads, dim(dads), loci.pos, chromo_mb, recom.rate)
         },
         stop("Invalid recombination type. Current options are 'map' or 'average'")
  )

  return(struct.rt)
}


#' Evolution of genetic and genomic landscapes
#'
#' This function simulates the evolution of individual-based populations forward in time.
#' @param x List of objects of class \code{"struct"} with the initial genetic structure of populations.
#' @param time Number of simulated generations.
#' @param type Type of simulated evolution. Current options are 'constant', 'dynamic', 'additive', and 'custom' (see Details).
#' @param recombination Type of recombination between homologous chromosomes. Current options are 'map' and 'average' (see Details).
#' @param recom.rate A numerical value for the recombination type 'average', or a vector of \eqn{nl - 1} elements, with the recombination rate between neighbouring loci for type 'map' (\eqn{nl}: number of loci) (see Details).
#' @param loci.pos A vector with the position of loci. The default value is NULL, but it is required for the recombination type 'average'.
#' @param chromo_mb A numerical value with the size of the simulated chromosome (in megabase). The default value is NULL, but it is required for the recombination type 'average'.
#' @param init.sex A list of vectors defining the sex of the initial individuals in the populations (1: females and 2: males). The default value is NULL and assigns the sex of the first generation randomly.
#' @param migration.rate A single value setting the migration rate between all populations, or a squared matrix of order equal to the number of populations, with the migration rate between them (migration['from', 'to']). It is ignored for single population simulations. The default  value is NULL with no migration between populations.
#' @param mutation.rate A numerical value setting the mutation rate per site. It is currently restricted to biallelic SNPs (genetic structures with values 1 or 2). The default value is NULL.
#' @param param.z A list with the parameter values for the function computing phenotypes. The list of parameters for each population should be included as a list of list (see Example)
#' @param param.w A list with the parameter values for the fitness function. The list of parameters for each population should be included as a list of list (see Example)
#' @param fun A character vector with the names of the custom phenotype and fitness functions. The default value is NULL, but it is required for the 'custom' type of evolution
#' @details This function returns a list of populations composed of two-dimensional arrays that represents a pair of homologous chromosomes. Rows represent individuals and columns the different loci. Each element of the array is an integer defining the copy of a given allele at a given locus.
#'
#' Different types of evolution are available for simulations:
#' \itemize{
#'  \item{'constant'} {a constant population size over time. There is no selection, equal sex ratio and each breeding pair generates two offspring. This case represent neutral evolution in which variation are due to recombination, mutation and migration. Variations in the population sizes are also possible due to the effect of migration.}
#'  \item{'dynamic'} {a dynamic population size over time. This type of evolution introduces to type 'constant' an unequal sex ratio, a variable number of offspring and a density-dependent demographic effect to avoid exponential growth. A list with the parameters for the phenotype function (param.z=list(sex.ratio)) and for the fitness function (param.w=list(mean.fitness, d.d)) are required:
#'    \itemize{
#'      \item{sex.ratio:} {a numerical value defining the sex ratio of populations. This value should be included within a list in 'param.z'.}
#'      \item{mean.fitness:} {an integer with the mean number of offsprings generated by breeding pair. This value should be included within a list in 'param.w'.}
#'      \item{d.d:} {numerical value introducing the density-dependent demographic effect to avoid exponential growth. This value should be included within a list in 'param.w'.}
#'
#'    The sex of individuals (1: females and 2: males) are generated by a binomial distribution with probability of being a male equal to the 'sex.ratio' The fitness of individuals (in number of offspring) are obtained from: \eqn{Poisson(\lambda) - N*d.d}, in which \eqn{\lambda} is equal to 'mean.fitness' and \eqn{N} represents the population size.
#'    }
#'   }
#'   \item{'additive'} {additive phenotype evolution. This function introduce quantitative phenotypes and convergent or divergent selection as implemented in Quilodrán et al. (2019). A list with the parameters for the phenotype function (param.z=list(sex.ratio, fitness.pos, bvs, add.loci, e.v)) and for the fitness function (param.w=list(b0,b1,b2,b3, d.v, add.loci)) are required.
#'
#'   The following parameters are needed for the computation of phenotypes (\eqn{z}):
#'    \itemize{
#'      \item{sex.ratio:} {A numerical value defining the sex ratio of populations. This value should be included within a list in 'param.z'.}
#'      \item{fitness.pos:} {A vector with the position of the additive loci participating in the computation of phenotypes. This value should be included within a list in 'param.z'.}
#'      \item{bvs:} {A matrix with the breeding values of alleles on each loci. The number of rows is equal to the number of additive loci, while the number of columns is equal to the maximum number of alleles in a locus. This object should be included within a list in 'param.z'.}
#'      \item{add.loci:} {An integer with the total number of additive loci participating in the computation of phenotypes. This object should be included within a list in 'param.z'.}
#'      \item{e.v:} {A numerical value defining a stochastic environmental variant in the computation of phenotypes. This value should be included within a list in 'param.z'.}
#'
#'      The default phenotype function (\eqn{z}) focus on an additive genetic genotype-phenotype map. Therefore, the sum of values of alleles at each locus gives a breeding value (\eqn{b_v}) for each individual at a given locus. The sum of breeding values \eqn{bvs} across loci gives a breeding value for the phenotypes (\eqn{z}), which is computed as follows:
#'      \deqn{z = \sum_{v =1 }^{n_a} b_v + \varepsilon_e(0, \sigma_e)}
#'      Where \eqn{n_a} is equal to 'add.loci' and \eqn{\sigma_e} is equal to 'e.v'. The environmental contribution \eqn{\varepsilon_e} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'e.v'.
#'      This function returns a data.frame with rows equal to the number of individuals and two columns ('sex' and 'z')
#'    }
#'   The default fitness function (\eqn{\omega}) computes a Gaussian relationship between \eqn{z} and \eqn{\omega}. The following parameters are needed:
#'    \itemize{
#'      \item{b0:} {A numerical value defining the maximum number of offspring generated by breeding pair. This value should be included within a list in 'param.w'.}
#'      \item{b1:} {A numerical value defining the phenotypic optima. In the gaussiam relationship between \eqn{z} and \eqn{\omega}, the value of 'b1' represent the \eqn{z} value expected to produce the maximum number of offspring. This value should be included within a list in 'param.w'. The difference in phenotypic optima between the populations drives the strength of 'divergent selection'. Populations exposed to equal phenotypic optima are considered to be under 'concordant selection'.}
#'      \item{b2:} {A numerical value defining the variance of the Gaussian curve. This value should be included within a list in 'param.z'.}
#'      \item{b3:} {A numerical value defining the intensity of the density-dependence on the fitness of individuals in a population of size \eqn{N}. This value should be included within a list in 'param.w'.}
#'      \item{d.v:} {A numerical value defining a stochastic demographic variant in the fitness of individuals. This value should be included within a list in 'param.w'.}
#'      \item{add.loci:} {An integer with the total number of additive loci participating in the computation of phenotypes. This object should be included within a list in 'param.w'.}
#'
#'      The default fitness function (\eqn{\omega}) has the form:
#'      \deqn{\omega = b_0  \exp^{ -\frac{1}{2} \left ( \frac{4z - b_1n_a}{b_2n_a} \right )^2} - b_3N + \varepsilon_d(0, \sigma_d )}
#'      Where \eqn{n_a} is equal to 'add.loci', \eqn{N} is the population size and \eqn{\sigma_d} is equal to 'd.v'. The demographic variant \eqn{\varepsilon_d} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'd.v'.
#'      This function returns a vector with the fitness value (\eqn{\omega}) of individuals.
#'    }
#'  }
#'  \item{'custom'} {is a custom computation of phenotypes and fitness functions.
#'  These functions will define the type of selection fitting particular case studies. The name of each user defined function should be introduced in the 'fun' parameter as a character vector with two elements e.g. c('phenotype', 'fitness').
#'  Assuming a per-generation time step, the potential number of offspring produced by each individual depends on its phenotype \eqn{\omega = f(z)}, which in turn depends on the individual genotype and on the environment \eqn{z = g(G, E)}.
#'  \eqn{G} is a numeric value determined by an individual’s genotype, representing the genetic value of the genotype. In the case of an additive genetic map, the genetic value of a genotype will be a breeding value. \eqn{E} represents the effect of the environment on phenotypic expression, and this enables the effects of plasticity on phenotypic expression to be captured.
#'  The list of parameters 'param.z' and 'param.w' include all needed parameters for the custom 'phenotype' and 'fitness' functions. These lists should not include variables.
#'  For the phenotype function, the variable genetic structure of individuals is already included as the object 'struct', which should also be the first argument of the custom phenotype function, followed by "...".
#'  This means that all parameters used in the custom phenotype function are only included in 'param.z'. The custom phenotype function should return a matrix with two columns, named "z" and "sex".
#'  The first column is the resulting phenotype value and the second column represents the assignation of sex to each individual. For the fitness function, the variable individual phenotype value, sex of individuals and the population size are already included in the environment.
#'  The function should start with these three elements ("z","sex","n"), followed by "...". The user does not need to use all three of these variables. All parameters needed for the custom fitness function should be included in 'param.w'.
#'  This function should return a vector 'w' with fitness values for the individuals (see Examples). }
#' }
#' The recombination between homologous chromosomes are either of type 'map' or 'average'. The first case needs a vector with the recombination rate (\eqn{\rho}) between neigbour loci of length equal to \eqn{nl - 1} (\eqn{nl}: number of loci). The probability of having a crossover (1) or not (0) is uniformly distributed at a rate defined by the value of \eqn{\rho} between loci
#' (i.e. positions with a probability smaller than \eqn{\rho} recombine). The uniform distribution allows each position with the same values of \eqn{\rho} to have an equal chance of crossover across all iterations. There is no recombination between homologous chromosomes when \eqn{\rho = 0}, both loci are completely linked (e.g. within an inversion or situated close to centromeres), while with a value of \eqn{\rho = 0.5}, the recombination rate is completely random (i.e. both loci are very distant on the same chromosome or are located on different chromosomes).
#' A value of \eqn{\rho < 0.5} means the loci are physically linked. The second case, when recombination is of type 'average', a numerical value with the average recombination rate per base pair should be supplied, the loci position and the size of the chromosome in megabase are also required. The crossover points are exponentially distributed as a Poisson process (see Example).
#' @return A list of objects of class \code{"struct"} or array.
#' @references
#' Quilodrán, C. S., Ruegg, K., Sendell-Price, A. T., Anderson, E., Coulson, T. and Clegg, S.  (2020).
#' The many population genetic and demographic routes to islands of genomic divergence.
#' \emph{Methods in Ecology and Evolution 11(1):6-21.}.
#' \doi{10.1111/2041-210X.13324}.
#' @seealso \code{\link{initial.struct}}
#' @examples
#' \dontrun{
#' ## We start with a population of 20 individuals and 10 biallelic SNPs
#' initial.population.size <- 20
#' n.loci <- 10
#' n.alleles.per.locus <- 2
#'
#' start1 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
#'
#' ##################################
#' # Type of evolution = "constant" #
#' ##################################
#'
#' ## We set a recombination map and the number of generations to be simulated
#' recom.map <- rep(0.5, n.loci-1) #all loci are independent
#' n.gens <- 10
#' pop <- evolve(list(start1), n.gens, "constant", "map", recom.map)
#'
#' ## A similar simulation but with recombination type "average".
#' ## We need to specify the position of loci and the size of the chromosome (MB)
#' loci.pos<- sample(80000: 12000000, 10) #random loci positions
#' chromo_mb<-12000000 	 #chromosome size
#' crossover <- 1/100000000.0 #average recombination rate (1 cM/MB)
#' pop <- evolve(list(start1), n.gens, "constant", "average", crossover, loci.pos, chromo_mb)
#'
#' # We include a mutation rate for the simulation of biallelic loci
#' pop <- evolve(list(start1), n.gens, "constant", "average", crossover,
#'               loci.pos, chromo_mb, mutation.rate = 0.0001)
#'
#' # A second population is included to incorporate the effect of migration
#' start2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
#' pop <- evolve(list(start1, start2), n.gens, "constant", "average", crossover,
#'               loci.pos, chromo_mb, migration.rate = 0.01)
#'
#' # The sex of individuals may be set for the starting generations.
#' sex.start1 <- sample(1:2, initial.population.size, replace=T)
#' sex.start2 <- sample(1:2, initial.population.size, replace=T)
#' init.sex <- list(sex.start1, sex.start2)
#' pop <- evolve(list(start1, start2), n.gens, "constant", "average", crossover,
#'               loci.pos, chromo_mb, init.sex = init.sex)
#'
#'
#' ##################################
#' # Type of evolution = "dynamic"  #
#' ##################################
#'
#' ## We set the parameters for the computation of phenotypes and for the
#' ## fitness function of the type 'dynamic'
#' sex.ratio <- 0.5
#' mean.fitness <- 3 #mean number of offsprings per breeding pair
#' d.d <- 0.01 #density-dependent demographic effect
#'
#' set.seed(1) #setting the seed for reproducible random numbers
#' pop <- evolve(list(start1), n.gens, "dynamic", "map", recom.map, param.z = list(list(sex.ratio)),
#'               param.w = list(list(mean.fitness, d.d)))
#'
#' ##################################
#' # Type of evolution = "additive" #
#' ##################################
#'
#' ## We set the parameters for the computation of phenotypes and for
#' ## the fitness function of the type 'additive'
#' sex.ratio <- 0.5
#' fitness.pos <- 1:n.loci #all simulated loci are additive
#' add.loci <- n.loci
#'
#' # next line is breeding value of additive loci
#' bvs <- t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus, add.loci)))
#'
#' e.v <- 0.01 					# stochastic environmental variant
#' param.z <- list(sex.ratio, fitness.pos, bvs, add.loci, e.v)
#'
#' b0 <- 6 # maximum number of offspring
#' b1 <- 1 # phenotypic optima
#' b2 <- 0.75 # variance of the Gaussian curve
#' b3 <- 0.01 # intensity of the density-dependence
#' d.v = 1 # stochastic demographic variant
#' param.w <- list(b0,b1,b2,b3, d.v, add.loci)
#'
#' set.seed(1) #setting the seed for reproducible random numbers
#' pop<-evolve(list(start1), n.gens, "additive", "map", recom.map,
#'             param.z = list(param.z), param.w = list(param.w))
#'
#' ##################################
#' #  Type of evolution = "custom"  #
#' ##################################
#'
#' ## We set a custom 'phenotype' and 'fitness' function.
#' ## In this example, the custom functions are very similar to the default 'additive' ones.
#'
#' phenotype2 <- function(struct, ...){
#'  pop.struct <- struct[ , fitness.pos, ]
#'  temp <- dim(pop.struct)
#'  mat <- matrix(1:temp[2],temp[1],temp[2],byrow=TRUE)
#'
#'  # the next line gives an array that gives the locus index at each position in pop.struct
#'  loci.n <- array(mat,c(temp[1],temp[2],2))
#'
#'  new.n <- (pop.struct-1)*add.loci+loci.n
#'  bvv <- as.vector(bvs) # turn to bvs
#'  outp <- array(bvv[new.n],c(temp[1],temp[2],2))
#'  z <- apply(outp,1,sum)-add.loci+rnorm(temp[1],0,e.v)
#'  sex <- rbinom(nrow(struct), 1, sex.ratio)+1
#'  return(cbind(z = z, sex = sex))
#' }
#'
#' fitness2 <- function(z, sex, n, ...){
#'  ng <- n.loci
#'  a <- b0
#'  b <- b1
#'  c <- b2
#'  w <- round( a*exp(-((z-b*ng)^2)/(2*(c*ng)^2) ) - b3*n + rnorm(n,0,d.v) , 0)
#'  w <- ifelse(w<0,0,w)
#'  return(w)
#' }
#'
#' set.seed(1) #setting the seed for reproducible random numbers
#' pop<-evolve(list(start1), n.gens, "custom", "map", recom.map,
#'            param.z =list(sex.ratio = sex.ratio, fitness.pos = fitness.pos,
#'                          bvs = bvs, add.loci = add.loci, e.v = e.v),
#'            param.w = list(b0 = b0, b1 = b1, b2 = b2, b3 = b3, d.v = d.v, n.loci = add.loci),
#'            fun=c("phenotype2", "fitness2"))
#'
#' }
#' @export
evolve <- function(x, time, type = c("constant", "dynamic", "additive", "custom"), recombination = c("map", "average"), recom.rate, loci.pos = NULL, chromo_mb = NULL, init.sex = NULL, migration.rate = NULL, mutation.rate = NULL, param.z=NULL, param.w=NULL, fun=c(phenotype=NULL, fitness=NULL)) {

  npop<-length(x)
  struct <- x

  if (is.list(x)==F) stop("Initial populations should be a list of arrays.", call.=F)
  if (recombination == "average" && (is.null(loci.pos) || is.null(chromo_mb))) stop("Loci position or size of the chromosome is missing with no default for 'average' recombination.", call.=F)
  if ((recombination == "map" && length(recom.rate) != (dim(x[[1]])[2] - 1)) || (recombination == "average" && length(recom.rate) != 1)) stop("Incorrect recombination rate for the type of recombination ('map' or 'average').", call.=F)

  if (type == "custom" && (is.character(fun) == FALSE || length(fun) != 2)) stop("One or both custom function names are missing or 'fun' is not well defined. 'fun' is a character vector of length 2, with the name of the phenotype and fitness functions e.g. c('phenotype', 'fitness').", call.=F)

   if (!is.null(mutation.rate)) {
    if (sum(sapply(1:npop, function(i) { length(table(x[[i]])) != 2 })) != 0) {
      mutation.rate = NULL
      warning("Mutation rate is ignored for non-biallelic genetic structures. The genetic structure input should contains integers of values 1 or 2.", call.=F)
    } else if (sum(sapply(1:npop, function(i) { 1 %in% x[[i]] && 2 %in% x[[i]] }) == FALSE) != 0) {
      mutation.rate = NULL
      warning("Mutation rate is ignored for non-biallelic genetic structures. The genetic structure input should contains integers of values 1 or 2.", call.=F)
    } }

  if (npop==1 && !is.null(migration.rate)) { migration.rate = NULL; warning("Migration is ignored for simulations with a single population.", call.=F) }
  if (!is.null(migration.rate)) {
    if (length(migration.rate) == 1) {
      disp=matrix(rep(migration.rate, npop*npop), nrow = npop, ncol = npop, byrow = TRUE)
      diag(disp)<-1
    } else {
      disp=migration.rate
      diag(disp)<-1
      if (sum(dim(migration.rate)!=npop) != 0) stop("Migration should be a single value or an square matrix of order equal to the number of populations.", call.=F)
    } }


  if (!is.null(init.sex) && (length(init.sex) != npop || sum(sapply(1: length(x), function(i){ nrow(x[[i]])!=length(init.sex[[i]]) }) ) != 0)) stop("There is not enough females or males assigned to the initial populations.", call.=F)

  ############################################
  pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total = time, clear = FALSE, width= 100)

  switch(recombination,
    map = {
      res <- list()
      for (i in 1:time) {
        if (i > 1) { init.sex <- NULL }

        y <- lapply(1:npop, function(i) { append(list(struct[[i]]), list(recom.rate, init.sex[[i]], mutation.rate, loci.pos = NULL, chromo_mb = NULL, param.z[[i]], param.w[[i]], fun)) })

        out <- lapply(y, newborns, recombination = recombination, type = type)

        if (!is.null(migration.rate)) {
          move <- combn(1:npop, 2)

          for (j in 1:ncol(move)) {
            n1 <- move[1, j]
            n2 <- move[2, j]
            rate1to2 <- disp[n1, n2]
            rate2to1 <- disp[n2, n1]
            outd<-migrate(out[[n1]], out[[n2]], rate1to2, rate2to1)
            out[[n1]]<-outd[[1]]
            out[[n2]]<-outd[[2]]
          }
        }

        struct<-out

        if (i %% time==0) res<- struct
        pb$tick()
      }
    },
    average = {
      res <- list()
      for (i in 1:time) {
        if (i > 1) { init.sex <- NULL }

        y <- lapply(1:npop, function(i) { append(list(struct[[i]]), list(recom.rate, init.sex[[i]], mutation.rate, loci.pos, chromo_mb, param.z[[i]], param.w[[i]], fun)) })

        out <- lapply(y, newborns, recombination = recombination, type = type)

        if (!is.null(migration.rate)) {
          move <- combn(1:npop, 2)

          for (j in 1:ncol(move)) {
            n1 <- move[1, j]
            n2 <- move[2, j]
            rate1to2 <- disp[n1, n2]
            rate2to1 <- disp[n2, n1]
            outd<-migrate(out[[n1]], out[[n2]], rate1to2, rate2to1)
            out[[n1]]<-outd[[1]]
            out[[n2]]<-outd[[2]]
          }
        }

        struct<-out

        if (i %% time==0) res<- struct
        pb$tick()
      }
    },
    stop("Invalid recombination type. Current options are 'map' or 'average'")
  )
  return(res)
}
