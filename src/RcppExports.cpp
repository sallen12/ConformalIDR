// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cidr_sequential
List cidr_sequential(NumericVector x_r, NumericVector x_out, List w, List W, List w_out, List pos_x, NumericVector y_unique_r, NumericVector y_out, int n_thr, int n_x);
RcppExport SEXP _ConformalIDR_cidr_sequential(SEXP x_rSEXP, SEXP x_outSEXP, SEXP wSEXP, SEXP WSEXP, SEXP w_outSEXP, SEXP pos_xSEXP, SEXP y_unique_rSEXP, SEXP y_outSEXP, SEXP n_thrSEXP, SEXP n_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_r(x_rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_out(x_outSEXP);
    Rcpp::traits::input_parameter< List >::type w(wSEXP);
    Rcpp::traits::input_parameter< List >::type W(WSEXP);
    Rcpp::traits::input_parameter< List >::type w_out(w_outSEXP);
    Rcpp::traits::input_parameter< List >::type pos_x(pos_xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_unique_r(y_unique_rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_out(y_outSEXP);
    Rcpp::traits::input_parameter< int >::type n_thr(n_thrSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    rcpp_result_gen = Rcpp::wrap(cidr_sequential(x_r, x_out, w, W, w_out, pos_x, y_unique_r, y_out, n_thr, n_x));
    return rcpp_result_gen;
END_RCPP
}
// cidr_static
List cidr_static(NumericVector x_r, NumericVector x_out, NumericVector w, List W, double w_out, List pos_x, NumericVector y_unique_r, NumericVector y_out, int n_thr, int n_x);
RcppExport SEXP _ConformalIDR_cidr_static(SEXP x_rSEXP, SEXP x_outSEXP, SEXP wSEXP, SEXP WSEXP, SEXP w_outSEXP, SEXP pos_xSEXP, SEXP y_unique_rSEXP, SEXP y_outSEXP, SEXP n_thrSEXP, SEXP n_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_r(x_rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_out(x_outSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< List >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type w_out(w_outSEXP);
    Rcpp::traits::input_parameter< List >::type pos_x(pos_xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_unique_r(y_unique_rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_out(y_outSEXP);
    Rcpp::traits::input_parameter< int >::type n_thr(n_thrSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    rcpp_result_gen = Rcpp::wrap(cidr_static(x_r, x_out, w, W, w_out, pos_x, y_unique_r, y_out, n_thr, n_x));
    return rcpp_result_gen;
END_RCPP
}
// lspm
arma::mat lspm(arma::colvec y_tr, arma::mat X_tr, arma::mat X_ts);
RcppExport SEXP _ConformalIDR_lspm(SEXP y_trSEXP, SEXP X_trSEXP, SEXP X_tsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type y_tr(y_trSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_tr(X_trSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_ts(X_tsSEXP);
    rcpp_result_gen = Rcpp::wrap(lspm(y_tr, X_tr, X_ts));
    return rcpp_result_gen;
END_RCPP
}
// olspm
arma::mat olspm(arma::colvec y_tr, arma::mat X_tr, arma::mat X_ts);
RcppExport SEXP _ConformalIDR_olspm(SEXP y_trSEXP, SEXP X_trSEXP, SEXP X_tsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type y_tr(y_trSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_tr(X_trSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_ts(X_tsSEXP);
    rcpp_result_gen = Rcpp::wrap(olspm(y_tr, X_tr, X_ts));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ConformalIDR_cidr_sequential", (DL_FUNC) &_ConformalIDR_cidr_sequential, 10},
    {"_ConformalIDR_cidr_static", (DL_FUNC) &_ConformalIDR_cidr_static, 10},
    {"_ConformalIDR_lspm", (DL_FUNC) &_ConformalIDR_lspm, 3},
    {"_ConformalIDR_olspm", (DL_FUNC) &_ConformalIDR_olspm, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ConformalIDR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
