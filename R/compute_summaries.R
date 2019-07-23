#' Conversion of class \code{struct} to class \code{loci}
#'
#' This function converts an object of class \code{"struct"} to an object of class \code{"loci"} that can be used by the package \code{pegas}
#' @param x List of objects of class \code{"struct"} with the initial genetic structure of populations.
#' @examples
#' ## We first create a random population with 100 individuals and 10 biallelic loci
#' N <- 100  # Population size
#' nl <- 10  # Number of additive loci
#' na <- 2  # Number of alleles per locus
#' G <- initial.struct(N,nl,na)
#'
#' ## We convert the object of class 'struct' into 'loci'
#' struct2pegas(list(G))
#' @export
struct2pegas <- function(x) {
  L <- dim(x[[1]])[2]
  npop <- length(x)

  pop<-lapply(1:npop, function(i){
    P <- paste(x[[i]][,,1], x[[i]][,,2], sep = "/") %>%
      matrix(ncol = L) %>%
      as.data.frame(stringsAsFactors = FALSE)
  })
  names(pop) <- 1:npop

  ret <- pop %>%
    dplyr::bind_rows(.id = "population")

  names(ret)[-1] <- paste("locus", 1:(ncol(ret)-1), sep = "_")

  # write it out so we can read it into pegas
  temp <- tempfile()
  write.table(ret, sep = " ", file = temp, col.names = TRUE, quote = FALSE, row.names = FALSE)

  pegas::read.loci(temp)
}
