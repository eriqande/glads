#include <Rcpp.h>
using namespace Rcpp;


//' Phenotype function from C++
//'
//' A function for the computation of additive phenotypes
//' @param G An object of class \code{"struct"} with the genetic structure of each individual in the population.
//' @param dims A vector with the dimensions of the 3-D array G.
//' @param bvs A matrix with the breeding values, indexed by loci and alleles
//' @param add_loci A vector with the position of additive loci participating in the computation of phenotypes.
//' @param sex_ratio A numerical value defining the sex ratio of populations
//' @param e_v A numerical value defining the environmental variance in the computation of phenotypes.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector rcpp_g2p_map(IntegerVector G, IntegerVector dims, NumericMatrix bvs, double add_loci, double sex_ratio, double e_v) {

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
    ret[i] = sum - add_loci;
  }

  ret = ret + rnorm(Ni, 0, e_v);

  retmat(_, 0) = ret;  // make the first column of the return matrix be z, the phenotype
  retmat(_, 1) = (runif(Ni) < sex_ratio) + 1;  // make the second column of the return be sex
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
//' @keywords internal
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
//' @keywords internal
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
//' crossing over points exponentially distributed as a Poisson process.
//' but that could be changed so that crossovers happen at a variable rate.   Note that this is hard-wired for diploidy.
//' @param G the structure giving the genotypes of the indviduals.  Actually a 3-D array indexed by indiv, locus, gene copy
//' @param dims the dimensions of the 3-D array G for internal use.
//' @param pos vector of positions of the loci.  This is an integer vector.  Has to be in sorted order (ascending)
//' @param chromo_length total chromoome length in base pairs
//' @param cross per base-pair rate of recombination.  For example, 1 cM per megabase equates to 1e-08.
//' @return  The return value is a long vector that can be squished into a matrix as appropriate to put it into
//' the genotype struct.
//' @export
//' @keywords internal
// [[Rcpp::export]]
IntegerVector rcpp_recombo_segregate_expo(IntegerVector G, IntegerVector dims, IntegerVector pos, int chromo_length, double cross) {
  int i,l,b, bl;  // for subscripting individual, locus
  int Ni = dims[0];
  int Nl = dims[1];
  int gam; // to say which gamete to choose (0 or 1)
  int nBreaks;
  int nb;


  IntegerVector breaks;  // to store recombination comparison variables
  IntegerVector ret(Ni * Nl, 0);  // allocate to the vector to return the gametes

  for(i=0; i<Ni; i++) {
    breaks = breakpoints1(chromo_length, cross);
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
//' @keywords internal
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

