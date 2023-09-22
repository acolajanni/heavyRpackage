#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//' Returns the inverse matrix using Rcpp Armadillo
//'
//' @export
// [[Rcpp::export]]
arma::mat invC(arma::mat M) {
  arma::mat A = inv(M) ;
  return A ;
 }


//'@rdname mvnpdf
//'@export
//'
// [[Rcpp::export]]
arma::vec mvnpdfoptimC(arma::mat x,
                        arma::colvec mean,
                        arma::mat varcovM,
                        bool Log=true){
  int p = x.n_rows;
  int n = x.n_cols;
  NumericVector y(n);

  mat Rinvsr = inv(trimatu(chol(varcovM)));
  double logSqrtDetvarcovM = sum(log(Rinvsr.diag()));
  double constant = - p * log2pi2;

  for (int i=0; i < n; i++) {
    colvec x_i = x.col(i) - mean;
    rowvec xRinvsr = trans(x_i) * Rinvsr;
    double quadform = sum(xRinvsr % xRinvsr);
    if (!Log) {
      y(i) = exp(-0.5*quadform + logSqrtDetvarcovM + constant);
    } else{
      y(i) = -0.5*quadform + logSqrtDetvarcovM + constant;
    }
  }

   return y;

 }

//'@rdname mvnpdf
//'@export
// [[Rcpp::export]]
arma::vec mvnpdfsmartC(arma::mat x,
                       arma::colvec mean,
                       arma::mat varcovM,
                       bool Log=true) {
  int p = x.n_rows;
  int n = x.n_cols;

  mat invvarcovM = inv_sympd(varcovM);
  double LogDetvarcovM = log_det_sympd(varcovM);

  vec y(n);
  for (int j=0; j < n; j++){
    colvec x_0 = x.col(j) - mean;
    vec tempres = -p/2.0 * log(2.0*M_PI) - 0.5 * LogDetvarcovM - 0.5 * x_0.t() * invvarcovM * x_0;
    y(j) =  tempres(0);
  }

  if (!Log){
    y =  exp(y);
  }

  return y;
}
