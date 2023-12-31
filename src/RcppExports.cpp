// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// invC
arma::mat invC(arma::mat M);
RcppExport SEXP _heavyRpackage_invC(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(invC(M));
    return rcpp_result_gen;
END_RCPP
}
// mvnpdfoptimC
arma::vec mvnpdfoptimC(arma::mat x, arma::colvec mean, arma::mat varcovM, bool Log);
RcppExport SEXP _heavyRpackage_mvnpdfoptimC(SEXP xSEXP, SEXP meanSEXP, SEXP varcovMSEXP, SEXP LogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type varcovM(varcovMSEXP);
    Rcpp::traits::input_parameter< bool >::type Log(LogSEXP);
    rcpp_result_gen = Rcpp::wrap(mvnpdfoptimC(x, mean, varcovM, Log));
    return rcpp_result_gen;
END_RCPP
}
// mvnpdfsmartC
arma::vec mvnpdfsmartC(arma::mat x, arma::colvec mean, arma::mat varcovM, bool Log);
RcppExport SEXP _heavyRpackage_mvnpdfsmartC(SEXP xSEXP, SEXP meanSEXP, SEXP varcovMSEXP, SEXP LogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type varcovM(varcovMSEXP);
    Rcpp::traits::input_parameter< bool >::type Log(LogSEXP);
    rcpp_result_gen = Rcpp::wrap(mvnpdfsmartC(x, mean, varcovM, Log));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _heavyRpackage_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_heavyRpackage_invC", (DL_FUNC) &_heavyRpackage_invC, 1},
    {"_heavyRpackage_mvnpdfoptimC", (DL_FUNC) &_heavyRpackage_mvnpdfoptimC, 4},
    {"_heavyRpackage_mvnpdfsmartC", (DL_FUNC) &_heavyRpackage_mvnpdfsmartC, 4},
    {"_heavyRpackage_timesTwo", (DL_FUNC) &_heavyRpackage_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_heavyRpackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
