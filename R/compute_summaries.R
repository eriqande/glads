
#' convert the two population "struct"s to a data frame of loci that the pegas package can use
#'
#' This will work as long as there are not more than 9 alleles.
#' @param P1 the pop struct, indexed by indiv, locus, gene copy, for pop 1
#' @param P2 the pop struct for pop 2
#' @export
struct2pegas <- function(P1, P2) {
  L <- dim(P1)[2]

  pop1 <- paste(P1[,,1], P1[,,2], sep = "/") %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)
  pop2 <- paste(P2[,,1], P2[,,2], sep = "/") %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)

  ret <- list("1" = pop1, "2" = pop2) %>%
    dplyr::bind_rows(.id = "population")

  names(ret)[-1] <- paste("locus", 1:(ncol(ret)-1), sep = "_")

  # write it out so we can read it into pegas
  temp <- tempfile()
  write.table(ret, sep = " ", file = temp, col.names = TRUE, quote = FALSE, row.names = FALSE)

  read.loci(temp)
}


#' convert the two population "struct"s to a data frame of loci that the pegas package can use
#'
#' This should work with any number of alleles
#' @param P1 the pop struct, indexed by indiv, locus, gene copy, for pop 1
#' @param P2 the pop struct for pop 2
#' @export
struct2pegas <- function(P1, P2) {
  L <- dim(P1)[2]

  pop1 <- paste(P1[,,1], P1[,,2], sep = "/") %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)
  pop2 <- paste(P2[,,1], P2[,,2], sep = "/") %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)

  ret <- list("1" = pop1, "2" = pop2) %>%
    dplyr::bind_rows(.id = "population")

  names(ret)[-1] <- paste("locus", 1:(ncol(ret)-1), sep = "_")

  # write it out so we can read it into pegas
  temp <- tempfile()
  write.table(ret, sep = " ", file = temp, col.names = TRUE, quote = FALSE, row.names = FALSE)

  pegas::read.loci(temp)
}



#' at each of the loci compute Fst
#'
#' uses pegas
#' @inheritParams struct2pegas
#' @export
fst_at_loci <- function(P1, P2) {
  dat <- struct2pegas(P1, P2)
  fst <- pegas::Fst(dat)
  dplyr::data_frame(locus = rownames(fst), Fst = fst[,"Fst"], Fit = fst[,"Fit"], Fis = fst[, "Fis"])

}


#' convert the two population "struct"s to a data frame of loci that the hierfstat package can use
#'
#' This will work as long as there are not more than 9 alleles. Probably will never end
#' up using this as hierfstat makes a very SLOW calculation of Fst.
#' @param P1 the pop struct, indexed by indiv, locus, gene copy, for pop 1
#' @param P2 the pop struct for pop 2
#' @export
struct2hierfstat <- function(P1, P2) {
  L <- dim(P1)[2]

  pop1 <- as.integer(paste(P1[,,1], P1[,,2], sep = "")) %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)
  pop2 <- as.integer(paste(P2[,,1], P2[,,2], sep = "")) %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)

  ret <- list("1" = pop1, "2" = pop2) %>%
    dplyr::bind_rows(.id = "population")

  names(ret)[-1] <- paste("locus", 1:(ncol(ret)-1), sep = "_")

  # write it out so we can read it into pegas
  temp <- tempfile()
  write.table(ret, sep = " ", file = temp, col.names = TRUE, quote = FALSE, row.names = FALSE)

  read.loci(temp)
}


