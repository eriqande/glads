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



//' simple function to return crossover points from an exponential
//'
//' we could swap this out for something more complicated with varying rates and
//' inversions fairly easily.
//' @param chr_len  the length of the chromosome in base pairs
//' @param rate  the rate of recombination per base pair (like 1/10^6)
//' @export
// [[Rcpp::export]]
IntegerVector breakpoints1(int chr_len, double rate) {

  int tot=0, expy;
  std::vector<int> ret;
  while(tot < chr_len) {
    expy = floor(rexp(1, rate)[0]);
    tot += expy;
    ret.push_back(tot);
  }
  return(wrap(ret));

}


//' rcpp version of function that does recombination and segregation with exponential crossovers
//'
//' In this version, we have to have a position for each locus (an integer less than 2^31) and then we have
//' crossing over points exponentially distributed as a Poisson process.  Mean crossover distance is 10^8 base pairs,
//' but that could be changed so that crossovers happen at a variable rate.   Note that this is hard-wired for diploidy.
//' @param G the structure giving the genotypes of the indviduals.  Actually a 3-D array indexed by indiv, locus, gene copy
//' @param dims the dimensions of the 3-D array G for internal use.
//' @param pos vector of positions of the loci.  This is an integer vector.  Has to be in sorted order (ascending)
//' @param chromo_length total chromoome length in base pairs
//' @return  The return value is a long vector that can be squished into a matrix as appropriate to put it into
//' the genotype struct.
//' @export
// [[Rcpp::export]]
IntegerVector rcpp_recombo_segregate_expo(IntegerVector G, IntegerVector dims, IntegerVector pos, int chromo_length) {
  int i,l,b, bl;  // for subscripting individual, locus
  int Ni = dims[0];
  int Nl = dims[1];
  int gam; // to say which gamete to choose (0 or 1)
  int nBreaks;
  int nb;


  IntegerVector breaks;  // to store recombination comparison variables
  IntegerVector ret(Ni * Nl, 0);  // allocate to the vector to return the gametes

  for(i=0; i<Ni; i++) {
    breaks = breakpoints1(chromo_length, 1.0/100000000.0);
    nBreaks = breaks.size(); /* store how many breaks there are */
    nb = nBreaks;

    gam = runif(1)[0] < 0.5;  // this sets the haplotype we start segregating from

    for(l=0; l<Nl; l++) {
      /* if number of breakoints that are > than the position changes, then we change gam. */
      for(b=0,bl=0 ;b<nBreaks; b++)  if(pos[l] < breaks[b]) bl++;
      if(bl < nb) {
        gam = (gam + (nb - bl)) % 2;  // note if two crossovers have occurred, there is no recomb.  nb-bl takes care of those even crossover numbers...
        nb = bl;
      }
      //Rcpp::Rcout << "indiv " << i << "  locus " << l << " pos[l] " << pos[l] << " gam " << gam << "  breaks  " << breaks << std::endl;
      ret[l*Ni + i] =  G[gam * (Ni*Nl) + l*(Ni) + i];
    }
  }

  return(ret);
}





//' dispersal function in rcpp
//'
//' a quick rcpp based implementation because the abind in the R implementation
//' gobbles up a lot of time.  For this you just pass in a vector for pop1 and a vector for
//' pop2 that says where each individual goes.  That way we can do the individual selection
//' outside of this function.
//' @param P1 first pop struct, indexed by indiv, locus, gene copy
//' @param P2 second pop struct
//' @param d1 dim of first pop struct
//' @param d2 dim of second pop struct
//' @param a1 assignments of individuals in pop 1 to either pop 1 or 2
//' @param a2 assignments of individuals in pop 2 to either pop 1 or 2
//' @export
// [[Rcpp::export]]
List rcpp_dispersal_placement(IntegerVector P1, IntegerVector P2, IntegerVector d1, IntegerVector d2, IntegerVector a1, IntegerVector a2) {

#define I(i,l,g,N)     g * (N * L) + l * N + i
  int i,l,g;
  int N1, N2;
  int r1 = 0, r2 = 0;  // for subscripting into the ret arrays

  int L = d1[1];  // number of loci
  int GC = 2;  // number of gene copies
  int I1 = d1[0]; // number of indivs in struct 1
  int I2 = d2[0];

  N1 = sum(a1 == 1) + sum(a2 == 1);
  N2 = sum(a1 == 2) + sum(a2 == 2);

  IntegerVector ret1(N1 * L * GC);
  IntegerVector ret2(N2 * L * GC);


  // shuffle off the critters from Pop1
  for(i=0;i<I1;i++) {
    if(a1[i]==1) {
      for(l=0;l<L;l++) {
        for(g=0;g<2;g++) ret1[I(r1,l,g,N1)] = P1[I(i,l,g,I1)];
      }
      r1++;
    } else {
      for(l=0;l<L;l++) {
        for(g=0;g<2;g++) ret2[I(r2,l,g,N2)] = P1[I(i,l,g,I1)];
      }
      r2++;
    }
  }

  // and then the ones from Pop2
  for(i=0;i<I2;i++) {
    if(a2[i]==1) {
      for(l=0;l<L;l++) {
        for(g=0;g<2;g++) ret1[I(r1,l,g,N1)] = P2[I(i,l,g,I2)];
      }
      r1++;
    } else {
      for(l=0;l<L;l++) {
        for(g=0;g<2;g++) ret2[I(r2,l,g,N2)] = P2[I(i,l,g,I2)];
      }
      r2++;
    }
  }

  return(List::create(ret1, ret2));
}

