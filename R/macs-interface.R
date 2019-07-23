
# functions to run the MaCS coalescent simulatar to initialize variation


#' path to the macs executable
#'
#' blab blah
#' @export
macs_path <- function() {
  if(.Platform$OS.type=="windows") {
    stop("Windows not supported")
  }
  else {
    if(.Platform$OS.type=="unix") {
      if(Sys.info()["sysname"]=="Darwin") {
        path <- file.path(system.file(package="glads"), "bin/macs")
      } else {
        stop("Linux not currently supported")
      }
    }
  }
  if(!file.exists(path)) stop(paste("macs executable should be installed at", path,"but does not seem to be there"))

  return(path)
}


#' path to the msformatter executable
#'
#' blab blah
#' @export
msformatter_path <- function() {
  if(.Platform$OS.type=="windows") {
    stop("Windows not supported")
  }
  else {
    if(.Platform$OS.type=="unix") {
      if(Sys.info()["sysname"]=="Darwin") {
        path <- file.path(system.file(package="glads"), "bin/msformatter")
      } else {
        stop("Linux not currently supported")
      }
    }
  }
  if(!file.exists(path)) stop(paste("msformatter executable should be installed at", path,"but does not seem to be there"))

  return(path)
}




#' prepare a macs command line to run with system
#'
#' this just writes that command line
#' @param num number of sequences to simulate
#' @param length total length of sequence in base pairs
#' @param A a list of arguments to pass to macs.  If you want just the flag, pass its value as "".  See the default for an idea
#' of what it should be.  Each element of the list should just be a string.
#' @param output name of the output file
#' @export
#' @examples
#' # a small short example.  Let's say we are going to simulate a sample of 200 sequences of
#' # a stretch of 10 million base
#' # pairs and we assume that the per-base-pair-per-generation neutral mutation rate is
#' # 1e-08, and that the recombination rate is roughly the same (1 cM per megabase).  Then our
#' # scaled rates (4Ne\mu and 4Ne\rho), assuming an effective size of 10,000 will be 4e-04.  This is
#' # what that looks like:
#' #
#' # first, the parameters
#' A <- list(t = 4e-04, r = 4e-04, h = 100)
#'
#' # then simulate it
#' system(macs_command_line(200, 1e7, A = A))
#'
#' # then read the results into a matrix
#' haps <- read_macs_output("haplotypes.txt")
#'
#' # check how many sequences and variants we have
#' dim(haps$mat)
#'
#' # see how many of those variants are singletons
#' sum(colSums(haps$mat-1) == 1)
macs_command_line <- function(num = 200,
                              length = 1e7,
                              A = list(
                                t = "0.0001",
                                r = "0.001",
                                h = "100"
                              ),
                              output = "haplotypes.txt")
{
  listy <- paste(paste("-", names(A), sep = ""), A, sep = " ", collapse = " ")
  call <- paste(macs_path(), num, length, listy, "2>trees_dump.txt |", msformatter_path(), ">", output)

  call
}


#' simulate starting variation for the simulation
#'
#' This is a high-level wrapper that calls \code{\link{macs_command_line}} and
#' \code{\link{read_macs_output}} and returns a bunch of biallelic genotypes that
#' one can then use however they want to initialize the forward in time simulation.
#' The parameterization here is built for convenience, for people that don't like
#' thinking about scaling things in coalescent time too much.  Basically you tell it
#' the length of the sequence in base pairs and the desired number of SNPs.  It assumes
#' by default that recombination rate is 1 cM per megabase.  You set the effective size
#' so that we can scale things as necessary by 4Ne.
#' @param n the number of sequences (2 times the number of diploid individuals you will want)
#' @param bp the number of base pairs in length of the segment to simulate
#' @param rho number of centiMorgans per megabase.  Default is 1.  (Note, we can change this
#' in the future to give ourselves recombination rate variation.)
#' @param Ne the coalescent effective size of the population over its coalescent history.
#' @param S the desired expected number of segregating sites. We will set theta to shoot for this.  Note that
#' doing so we can have high recombination across an entire chromosome, but simulate getting SNPs from
#' only a fraction of the total area of the chromosome.
#' @param h how many base pairs back to extend the Markovian approximation using macs
#' @export
coalescent_starting_variation <- function(n = 100, bp = 1e6, rho = 1, Ne = 1e4, S = 1000, h = 100) {
  if(FALSE) {
    n <- 100
    bp <- 1e6
    rho <- 1
    Ne <- 1e4
    S <- 1000
    h <- 100
  }

  # figure out what we would like theta across the whole segment to be.
  hm <- sum(1/1:(n-1))
  theta <- S/hm
  t <- theta / bp  # this is the per-base-pair t argument to macs

  # figure out the per base pair recombination rate
  r <- rho * 10^-8 * Ne

  # make an args list out of that
  A <- list(t = t, r = r, h = h)

  # do the call
  system(macs_command_line(n, bp, A = A))

  # slurp up the results and change the positions to base pairs
  ret <- read_macs_output("haplotypes.txt")

  ret$pos <- floor(ret$pos * bp)

  ret
}
