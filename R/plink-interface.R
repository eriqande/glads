

# functions for doing stuff from within plink. Mostly this will be
# writing input files, computing Fst, etc.

#' write plink .ped and map files from simulated markers
#'
#' The just assumes that everything is on chromosom 1.
#' @param R a list of struct1 and struct2---the structures holding the alleles
#' carried by the indivs
#' @param pos a vector giving the positions (in base pairs) of the loci
#' @param plinkout the prefix of the file to write out.
result2plink <- function(R, pos, plinkout = "plink") {
  # get the genotype matrices first
  x1 <- R[[1]]
  y1 <- array(paste(x1[,,1], x1[,,2]), dim = dim(x1)[-length(dim(x1))])

  x2 <- R[[2]]
  y2 <- array(paste(x2[,,1], x2[,,2]), dim = dim(x2)[-length(dim(x2))])

  # then add the family and indiv ID columns and such for plink format.
  # We will let population 1 be "Unaffecteds" and population 2 be "Affecteds"
  # so we can easily say how to compute Fst without using a "within" file.
  z1 <- cbind(paste("a", 1:nrow(y1), sep=""),
              paste("a", 1:nrow(y1), sep=""),
              0,
              0,
              0,
              1,
              y1)
  z2 <- cbind(paste("b", 1:nrow(y2), sep=""),
              paste("b", 1:nrow(y2), sep=""),
              0,
              0,
              0,
              2,
              y2)

  write.table(rbind(z1, z2),
              file = paste(plinkout, ".ped", sep = ""),
              sep = "  ",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  # now, write out the map file:
  map <- cbind(1,
               paste("snp", 1:length(pos), sep = ""),
               0,
               pos)

  write.table(map,
              file = paste(plinkout, ".map", sep = ""),
              sep = "  ",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}



#' return a data frame of Fst values given a struct
#'
#' blah blah
#' @inheritParams result2plink
#' @export
fst_via_plink <- function(R, pos) {
  tmp <- "xx-plink-fst"
  # first, remove any other plink output files that might be there
  unlink(paste(tmp, c("ped", "map", "fst"), sep = "."))

  # then write new ones
  result2plink(R, pos, plinkout = tmp)

  # then run plink and compute Fst
  system("plink -file xx-plink-fst --fst case-control --allow-no-sex --out xx-plink-fst > plink_dump")

  # then read the output into a tibble (and return that!)
  suppressMessages(readr::read_delim("xx-plink-fst.fst", delim = "\t"))
}
