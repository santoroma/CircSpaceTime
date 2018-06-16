// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  // Thanks to Ahmadou Dicko http://gallery.rcpp.org/articles/simulate-multivariate-normal/
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
