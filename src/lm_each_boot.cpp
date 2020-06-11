#include <Rcpp.h>
using namespace Rcpp;

//' lm_each_boot
//' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
//' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
//' @param n number of multiprobability
// [[Rcpp::export]]

Rcpp::List lm_each_boot1(Rcpp::Formula formula, Rcpp::DataFrame data, int n){
  Rcpp::Environment myEnv = Rcpp::Environment::global_env();
  Rcpp::Environment base("package:stats");
  Function  lm1= myEnv["lm1"];
  Rcpp::Function lm_rmu = base["rmultinom"];
  int a = data.nrows();
  NumericVector freqs = lm_rmu(1,n,rep(1,a));
  Rcpp::List lm2=lm1(formula,data,freqs);
  return lm2;
}



