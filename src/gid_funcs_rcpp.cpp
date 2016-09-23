#include <Rcpp.h>
using namespace Rcpp;



//' rcpp version of g2p.map to see how it fares speedwise
//'
//' Just a reimplementation using Rcpp
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

  for(i=0; i<Ni; i++) {
    sum = 0.0;
    for(l=0; l<Nl; l++) {
      for(a=0; a<Na; a++) {
        sum += bvs(l, G[a*(Ni*Nl) + l*(Ni) + i] - 1);
      }
    }
    ret[i] = sum - num_loci_t;
  }

  return(ret + rnorm(Ni, 0, v_e));
}
