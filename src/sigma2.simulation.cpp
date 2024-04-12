#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title Simulate sigma squared
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters of the model
//' @param data Dataset of type matrix
//' @param parm Current estimated parameters
//' @return Sample of idiosyncratic variances, sigma squared
//' @noRd
// [[Rcpp::export]]


arma::vec sigma2_simulation(List attributes, List hyp_parm, arma::mat data, List parm){
  
  int R = attributes["R"];
  int N = attributes["N"];
  double n_sigma = hyp_parm["n.sigma"];
  double n_s2_sigma = hyp_parm["n.s2.sigma"];
  mat B = parm["B"];
  mat X = parm["X"];
  vec sigma2out(R);
  
  for(int j = 0; j < R; j++){
    double n_sigma2_aux = n_sigma + N;
    double n_s2_sigma2_aux = n_s2_sigma;
    for(int i = 0; i < N; i++){
      n_s2_sigma2_aux += pow(data(i,j) - conv_to<double>::from(B.row(j) * conv_to<colvec>::from(X.row(i))), 2.0);
    }
    sigma2out(j) = 1 / as<double>(wrap(R::rgamma(n_sigma2_aux/2, 2/n_s2_sigma2_aux)));
  }
  
  return sigma2out;
  
}

