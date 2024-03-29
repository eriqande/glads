% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glads_funcs.R
\name{fitness}
\alias{fitness}
\title{Fitness function}
\usage{
fitness(z, N, b0, b1, b2, b3, d.v, add.loci)
}
\arguments{
\item{z}{A data.frame with the phenotype value ('z') and the sex of each individual in a population.}

\item{N}{Population size of the population.}

\item{b0}{A numerical value defining the maximum number of offspring that may be generated by a breeding pair.}

\item{b1}{A numerical value defining the phenotypic optima of a given population (see Details).}

\item{b2}{A numerical value defining the variance of the Gaussian curve.}

\item{b3}{A numerical value defining the intensity of the density-dependence on the fitness of individuals in a population of size 'N'.}

\item{d.v}{A numerical value defining a stochastic demographic variant in the fitness of individuals.}

\item{add.loci}{An integer with the total number of additive loci participating in the computation of phenotypes.}
}
\value{
A vector with the fitness value (\eqn{\omega}) for each individual.
}
\description{
This function computes a fitness value (\eqn{\omega}) depending on the phenotype (z) of individuals. The relationship between both variables is assumed to be Gaussian.
}
\details{
This function is used internally in the function evolve() of type 'selection' in order to compute the fitness value of individuals.

The value of 'b1' represents the \eqn{z} value expected to produce the maximum number of offspring.
The difference in phenotypic optima between the populations drives the strength of 'divergent selection'. Populations exposed to equal phenotypic optima are considered to be under 'concordant selection'.

This fitness function (\eqn{\omega}) has the following form:
\deqn{\omega = b_0  \exp^{ -\frac{1}{2} \left ( \frac{4z - b_1n_a}{b_2n_a} \right )^2} - b_3N + \varepsilon_d(0, \sigma_d )}
Where \eqn{n_a} is equal to 'add.loci', \eqn{N} is the population size and \eqn{\sigma_d} is equal to 'd.v'. The demographic variant \eqn{\varepsilon_d} is assumed to be stochastic and normally distributed, with a mean of 0 and standard variation 'd.v'.
}
\examples{
## We first create a random population with 100 individuals and 10 loci
N <- 100  # Population size
nl <- 10  # Number of additive loci
na <- 4  # Number of alleles per locus
G <- initial.struct(N,nl,na)

## Additional parameters are needed for the computation of phenotypes
bvs <- t(array( seq(0,1, length = na) ,c(na, nl)))
sex.ratio <- 0.5
e.v=0.01

## Now we compute the additive phenotype value of individuals
phen <- phenotype(G, bvs, nl, sex.ratio, e.v)

## Additional parameters are needed for the fitness function
b0 <- 6
b1 <- 0.25
b2 <- 0.5
b3 <- 0.01
d.v = 1

## The fitness values of individuals are computed as follows:
fitness(phen, N, b0, b1, b2, b3, d.v, nl)
}
\references{
Quilodrán, C. S., Ruegg, K., Sendell-Price, A. T., Anderson, E., Coulson, T. and Clegg, S.  (2020).
The multiple population genetic and demographic routes to islands of genomic divergence.
\emph{Methods in Ecology and Evolution}.
\doi{10.1111/2041-210X.13324}.
}
\seealso{
\code{\link{evolve}} \code{\link{phenotype}}
}
