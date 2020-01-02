#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

double f(double t) {
  return exp(-abs(t));
}

//' @title Metropolis sampling
//' @description Use cpp to implement Metropolis sampling
//' @import Rcpp
//' @importFrom Rcpp evalCpp
//' @importFrom stats rnorm
//' @examples
//' dontrun{
//'   x0 = 25
//'   N = 2000
//'   sigma = 2
//'   rwC=rwMetropolis(sigma,x0,N)
//'   rwc
//' }
//' @export
// [[Rcpp::export]]
NumericVector rMetropolis (double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0;
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (f(y[0]) / f(x[i-1]))){
      x[i] = y[0];
    }
    else {
      x[i] = x[i-1];
    }
  }
  return(x);
}

