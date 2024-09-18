// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fit_lspm(arma::colvec y_tr, arma::mat X_tr, arma::mat X_ts) {
  arma::mat XX = inv(X_tr.t() * X_tr);
  arma::mat H = X_tr * XX * X_tr.t();

  int n = X_tr.n_rows;
  int m = X_ts.n_rows;
  arma::mat C(m, n);
  for (int j = 0; j < m; j++) {
    arma::colvec g = join_cols(X_tr * XX * X_ts.row(j).t(), X_ts.row(j) * XX * X_ts.row(j).t());
    arma::mat Hbar(n + 1, n + 1);
    Hbar.submat(0, 0, n - 1, n - 1) = H - g.head(n) * g.head(n).t() / (1 + g(n));
    Hbar.col(n) = g / (1 + g(n));
    Hbar.row(n) = Hbar.col(n).t();
    arma::colvec B = Hbar.submat(0, 0, n - 1, n).col(n) / sqrt(1 - Hbar.submat(0, 0, n - 1, n - 1).diag());
    B += sqrt(1 - Hbar(n, n));
    arma::colvec A = repmat((Hbar.submat(0, 0, n - 1, n).col(n).t() * y_tr) / sqrt(1 - Hbar(n, n)), n, 1);
    A += (y_tr - sum(Hbar.submat(0, 0, n - 1, n - 1) % repmat(y_tr.t(), n, 1), 1)) / sqrt(1 - Hbar.submat(0, 0, n - 1, n - 1).diag());
    C.row(j) = sort((A / B).t());
  }

  return C;
}


// [[Rcpp::export]]
arma::mat fit_olspm(arma::colvec y_tr, arma::mat X_tr, arma::mat X_ts) {
  arma::mat XX = inv(X_tr.t() * X_tr);
  arma::mat H = X_tr * XX * X_tr.t();

  int n = X_tr.n_rows;
  int m = X_ts.n_rows;
  arma::mat C(m, n);
  for (int j = 0; j < m; j++) {
    arma::colvec g = join_cols(X_tr * XX * X_ts.row(j).t(), X_ts.row(j) * XX * X_ts.row(j).t());
    arma::mat Hbar(n + 1, n + 1);
    Hbar.submat(0, 0, n - 1, n - 1) = H - g.head(n) * g.head(n).t() / (1 + g(n));
    Hbar.col(n) = g / (1 + g(n));
    Hbar.row(n) = Hbar.col(n).t();

    arma::colvec B = 1 - Hbar(n, n) + Hbar.submat(0, 0, n - 1, n).col(n);
    arma::colvec A = repmat((Hbar.submat(0, 0, n - 1, n).col(n).t() * y_tr), n, 1);
    A += (y_tr - sum(Hbar.submat(0, 0, n - 1, n - 1) % repmat(y_tr.t(), n, 1), 1));
    C.row(j) = sort((A / B).t());
  }

  return C;
}
