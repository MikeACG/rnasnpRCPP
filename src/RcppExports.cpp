// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_bpd
List compute_bpd(NumericMatrix wild, NumericMatrix mut, int regionX, int regionY);
RcppExport SEXP _rnasnpRCPP_compute_bpd(SEXP wildSEXP, SEXP mutSEXP, SEXP regionXSEXP, SEXP regionYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type wild(wildSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mut(mutSEXP);
    Rcpp::traits::input_parameter< int >::type regionX(regionXSEXP);
    Rcpp::traits::input_parameter< int >::type regionY(regionYSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_bpd(wild, mut, regionX, regionY));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _rnasnpRCPP_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rnasnpRCPP_compute_bpd", (DL_FUNC) &_rnasnpRCPP_compute_bpd, 4},
    {"_rnasnpRCPP_rcpp_hello_world", (DL_FUNC) &_rnasnpRCPP_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_rnasnpRCPP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}