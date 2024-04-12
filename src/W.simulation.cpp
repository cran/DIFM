#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title Simulate Ws
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters of the model
//' @param parm Current estimated parameters
//' @return Sample of evolution variance, W
//' @noRd
// [[Rcpp::export]]

arma::mat W_simulation(List attributes, List hyp_parm, List parm){
  
  // int L = attributes["L"];
  int N = attributes["N"];
  double n_w = hyp_parm["n.w"];
  mat GG = attributes["GG"];
  mat theta = parm["theta"];
  int LL = GG.n_cols;
  mat Psi = hyp_parm["Psi"];
  mat Wout(LL, LL);
  
  double n_w_aux = n_w + N - 1;
  mat n_s2_w_aux(LL, LL);
  for(int i = 1; i < N; i++){
    colvec thetanow = conv_to<colvec>::from(theta.row(i)) - GG * conv_to<colvec>::from(theta.row(i-1));
    n_s2_w_aux += thetanow * trans(thetanow);
  }
  Wout = arma::iwishrnd(Psi + n_s2_w_aux, n_w_aux);
  
  return Wout;

  
}
