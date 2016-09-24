library(microbenchmark)

# The code here runs a simulation to explore how a range of processes can generate genomic patterns.
# Motivated by the ideas in Sonya and Kristen's NERC and John Fell Fund grants.

# This chunk of code will likely be moved to a different file, with the following functions that do the hard
# work put in a separate file that is called from this one.

# First, clear workspace and graphics
rm(list=ls())
graphics.off()

# Higher level parameters that the user will set to run the simulation
initial.population.size <- 50 # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- 100 # how many linked genes we will deal with
n.alleles.per.locus     <- 5   # we will assume that each locus has the same number of alleles to start with
bv.for.alleles          <- t(array(1:5,c(n.alleles.per.locus,n.loci)))/4-0.25
V.e                     <- 1 # standard deviation in environmental component of the phenotype


# fire it up.
struct <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)
zs         <- g2p.map(struct,bv.for.alleles,n.loci,V.e)
mean(zs)



# benchmark the Rcpp version of g2p_map.  About a 30x speedup
mbres <- microbenchmark(zs  <- g2p.map(struct,bv.for.alleles,n.loci,V.e),
               z2 <- rcpp_g2p_map(struct, dim(struct), bv.for.alleles, n.loci, V.e),
              times = 1000)



# ensure that we get the same result with each function:
set.seed(1)
r_version <- lapply(1:100, function(x) g2p.map(struct,bv.for.alleles,n.loci,V.e))
set.seed(1)
rcpp_version <- lapply(1:100, function(x) rcpp_g2p_map(struct, dim(struct), bv.for.alleles, n.loci, V.e))
all.equal(r_version, rcpp_version)
