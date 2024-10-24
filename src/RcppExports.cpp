// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// corC
double corC(NumericVector x, NumericVector y);
RcppExport SEXP _fast_EOT_corC(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(corC(x, y));
    return rcpp_result_gen;
END_RCPP
}
// fastCor
NumericVector fastCor(NumericMatrix x_i, NumericVector y_i, bool standardised);
RcppExport SEXP _fast_EOT_fastCor(SEXP x_iSEXP, SEXP y_iSEXP, SEXP standardisedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_i(x_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_i(y_iSEXP);
    Rcpp::traits::input_parameter< bool >::type standardised(standardisedSEXP);
    rcpp_result_gen = Rcpp::wrap(fastCor(x_i, y_i, standardised));
    return rcpp_result_gen;
END_RCPP
}
// lmC
List lmC(NumericVector x, NumericVector y);
RcppExport SEXP _fast_EOT_lmC(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(lmC(x, y));
    return rcpp_result_gen;
END_RCPP
}
// respLmParam
List respLmParam(NumericMatrix x, NumericMatrix y, int cell);
RcppExport SEXP _fast_EOT_respLmParam(SEXP xSEXP, SEXP ySEXP, SEXP cellSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type cell(cellSEXP);
    rcpp_result_gen = Rcpp::wrap(respLmParam(x, y, cell));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fast_EOT_corC", (DL_FUNC) &_fast_EOT_corC, 2},
    {"_fast_EOT_fastCor", (DL_FUNC) &_fast_EOT_fastCor, 3},
    {"_fast_EOT_lmC", (DL_FUNC) &_fast_EOT_lmC, 2},
    {"_fast_EOT_respLmParam", (DL_FUNC) &_fast_EOT_respLmParam, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fast_EOT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
