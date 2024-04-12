#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title Simulate Xs
//' @param attributes Attributes of the model and the data
//' @param parm Current estimated parameters
//' @return Sample of common factors, X
//' @noRd
// [[Rcpp::export]]

arma::mat X_simulation(List attributes, List parm){
  
  int N = attributes["N"];
  int L = attributes["L"];
  mat FF = attributes["FF"];
  mat theta = parm["theta"];
  mat Xout(N,L);
  
  for(int i = 0; i < N; i++){
    Xout.row(i) = conv_to<rowvec>::from(FF * conv_to<colvec>::from(theta.row(i)));
  }
  
  return Xout;
  
}
