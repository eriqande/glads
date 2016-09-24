#include <Rcpp.h>
using namespace Rcpp;



//' rcpp version of g2p.map to see how it fares speedwise
//'
//' Just a reimplementation using Rcpp.  This returns a matrix.  The first column
//' are the phenotype values and the second column are the sexes.
//' @param G the structure giving the genotypes of the indviduals.  Actually a 3-D array indexed by indiv, locus, gene copy
//' @param dims the dimensions of the 3-D array G for internal use.
//' @param bvs matrix of the breeding values, indexed by loci and alleles
//' @param num_loci_t  the number of loci, essentially something to get the mean Z to be 0
//' @param v_e the environmental variance in the phenotype
//' @export
// [[Rcpp::export]]
NumericVector rcpp_g2p_map(IntegerVector G, IntegerVector dims, NumericMatrix bvs, double num_loci_t, double v_e) {

  int i,l,a;  // for subscripting individual, locus, gene copy
  int Ni = dims[0];
  int Nl = dims[1];
  int Na = dims[2];
  double sum;

  // allocate space to return the result
  NumericVector ret(Ni, 0.0);
  NumericMatrix retmat(Ni, 2);

  for(i=0; i<Ni; i++) {
    sum = 0.0;
    for(l=0; l<Nl; l++) {
      for(a=0; a<Na; a++) {
        sum += bvs(l, G[a*(Ni*Nl) + l*(Ni) + i] - 1);
      }
    }
    ret[i] = sum - num_loci_t;
  }

  ret = ret + rnorm(Ni, 0, v_e);

  retmat(_, 0) = ret;  // make the first column of the return matrix be z, the phenotype
  retmat(_, 1) = (runif(Ni) < 0.5) + 1;  // make the second column of the return be sex
  colnames(retmat) = CharacterVector::create("z", "sex");

  return(retmat);
}



//' rcpp version of function that does recombination and segregation
//'
//' Note that this is hard-wired for diploidy.
//' @param G the structure giving the genotypes of the indviduals.  Actually a 3-D array indexed by indiv, locus, gene copy
//' @param dims the dimensions of the 3-D array G for internal use.
//' @param rf vector of recombination fractions.  There should be one minus the number of loci. The first
//' one corresponds to recombination between the first and the second marker.  These are the probabilities
//' of a recombination between the markers during a meiosis.
//' @return  The return value is a long vector that can be squished into a matrix as appropriate to put it into
//' the genotype struct.
//' @export
// [[Rcpp::export]]
IntegerVector rcpp_recombo_segregate(IntegerVector G, IntegerVector dims, NumericVector rf) {
  int i,l;  // for subscripting individual, locus
  int Ni = dims[0];
  int Nl = dims[1];
  int gam; // to say which gamete to choose (0 or 1)

  NumericVector rando(Nl);  // to store recombination comparison variables
  IntegerVector ret(Ni * Nl, 0);  // allocate to the vector to return the gametes

  for(i=0; i<Ni; i++) {
    rando = runif(Nl);
    gam = rando[0] < 0.5;  // this sets the haplotype we start segregating from
    for(l=0; l<Nl; l++) {
      if(l > 0) { // if we are not on the first locus, allow the chance for a recombination
        gam = (gam + (rando[l] < rf[l-1])) % 2; // if recombination, change gam from 0 to 1, or 1 to 0.
      }
      //Rcpp::Rcout << "indiv " << i << "  locus " << l << " rando[l] " << rando[l] << "  rf[l-1] " << rf[l-1] << "  gam: " << gam << std::endl;
      ret[l*Ni + i] =  G[gam * (Ni*Nl) + l*(Ni) + i];
    }
  }

  return(ret);
}
