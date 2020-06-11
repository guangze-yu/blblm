#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List lm_each_boot1(Rcpp::Formula formula, Rcpp::DataFrame data, int n){
  Rcpp::Environment myEnv = Rcpp::Environment::global_env();
  Rcpp::Environment base("package:stats");
  Function  lm1= myEnv["lm1"];
  Rcpp::Function lm_rmu = base["rmultinom"];
  int a = data.nrows();
  Rcpp::List freqs = lm_rmu(1,n,rep(1,a));
  return lm1(formula,data,freqs);
}



