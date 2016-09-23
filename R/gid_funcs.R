#' a function to generate the initial population structure
#'
#' fill in initial values
#' @param n.N1 Number of individuals in population 1
#' @param n.l Number of loci
#' @param n.a.l Number of alleles at each locus.  Just set as a single scalar for all the loci right now.
#' @export
initial.struct <- function(n.N1, n.l, n.a.l) {
  rand.ints           <- sample(1:n.a.l, n.N1 * n.l * 2, replace = TRUE) # alleles at each locus equifrequent
  struct              <- array(rand.ints, c(n.N1, n.l, 2)) # a 3-D array
  return(struct)
}

#' make the phenotype
#'
#' This maps genotype to phenotype
#' @param pop.struct a structure holding the population members.  A 3-D array subscripted by indivs, loci, gene-copies
#' @param bvs array of breeding values.  an array with n.loci rows and number of columns equal to number of alleles.
#' @param n.loci.t
#' @export
g2p.map <- function(pop.struct, bvs, n.loci.t, ve){
    temp <- dim(pop.struct)
    mat <- matrix(1:temp[2],temp[1],temp[2],byrow=TRUE)
    loci.n <- array(mat,c(temp[1],temp[2],2)) # array that gives the locus index at each position in pop.struct
    new.n <- (pop.struct-1)*n.loci.t+loci.n
    bvv <- as.vector(bvs) # turn to bvs
    outp <- array(bvv[new.n],c(temp[1],temp[2],2))
    z <- apply(outp,1,sum)-n.loci.t+rnorm(temp[1],0,ve)
    return(z)
}

#' estimate fitness given phenotype
#' @export
fitness <- function(z,n,sigma){
      w <- round(0.3 + 0.1*z -0.002*n + rnorm(n,0,0.5))
}

#' makes a gamete from a genotype
#'
#' the genotype is given as a matrix of alleles.  First column is the first haplotype and the second
#' column is the second haplotype.
#' @export
recombine <- function(genotype){
  temp <- dim(genotype)
  i <- sample(c(1,2),1,FALSE)
  j <- if (i==1) j <- 2 else j <- 1
  x <- c(genotype[,i],genotype[,j])
  ra <- runif(temp[1]) # this is the probability of recombining at a locus
  re <- ifelse(ra>0.95,1,0)
  re <- cumsum(re)
  re <- (re %% 2)+1
  xx <- 1:temp[1]
  gamete <- x[ifelse(re==1,xx,xx+100)]
  return(gamete)
}

