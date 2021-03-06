% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glads_funcs.R
\name{mutation}
\alias{mutation}
\title{Mutation function}
\usage{
mutation(G, mutation.rate)
}
\arguments{
\item{G}{An object of class \code{"struct"} with the genetic structure of individuals before gametogenesis.}

\item{mutation.rate}{Mutation rate of biallelic loci}
}
\value{
An object of class \code{"struct"} or an array.
}
\description{
This function introduces mutations on biallelic genetic structures.
}
\details{
This function is used internally in the function newborns() in order to introduce mutations on biallelic genetic structures before meiosis. The genotypes are introduced as an array of integers defining the copy of a given allele at a given locus. Rows represent individuals and columns the different loci.
}
\seealso{
\code{\link{evolve}} \code{\link{newborns}}
}
\keyword{internal}
