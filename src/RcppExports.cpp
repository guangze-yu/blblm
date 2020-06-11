// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lm_each_boot1
Rcpp::List lm_each_boot1(Rcpp::Formula formula, Rcpp::DataFrame data, int n);
RcppExport SEXP _blblm_lm_each_boot1(SEXP formulaSEXP, SEXP dataSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Formula >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(lm_each_boot1(formula, data, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_blblm_lm_each_boot1", (DL_FUNC) &_blblm_lm_each_boot1, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_blblm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
