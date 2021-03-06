% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_recombo_segregate_expo}
\alias{rcpp_recombo_segregate_expo}
\title{rcpp version of function that does recombination and segregation with exponential crossovers}
\usage{
rcpp_recombo_segregate_expo(G, dims, pos, chromo_length, cross)
}
\arguments{
\item{G}{the structure giving the genotypes of the indviduals.  Actually a 3-D array indexed by indiv, locus, gene copy}

\item{dims}{the dimensions of the 3-D array G for internal use.}

\item{pos}{vector of positions of the loci.  This is an integer vector.  Has to be in sorted order (ascending)}

\item{chromo_length}{total chromoome length in base pairs}

\item{cross}{per base-pair rate of recombination.  For example, 1 cM per megabase equates to 1e-08.}
}
\value{
The return value is a long vector that can be squished into a matrix as appropriate to put it into
the genotype struct.
}
\description{
In this version, we have to have a position for each locus (an integer less than 2^31) and then we have
crossing over points exponentially distributed as a Poisson process.
but that could be changed so that crossovers happen at a variable rate.   Note that this is hard-wired for diploidy.
}
\keyword{internal}
