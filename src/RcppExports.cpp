// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_g2p_map
NumericVector rcpp_g2p_map(IntegerVector G, IntegerVector dims, NumericMatrix bvs, double num_loci_t, double v_e);
RcppExport SEXP gids_rcpp_g2p_map(SEXP GSEXP, SEXP dimsSEXP, SEXP bvsSEXP, SEXP num_loci_tSEXP, SEXP v_eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bvs(bvsSEXP);
    Rcpp::traits::input_parameter< double >::type num_loci_t(num_loci_tSEXP);
    Rcpp::traits::input_parameter< double >::type v_e(v_eSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_g2p_map(G, dims, bvs, num_loci_t, v_e));
    return rcpp_result_gen;
END_RCPP
}