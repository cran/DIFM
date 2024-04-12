#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title Simulate taus
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters of the model
//' @param parm Current estimated parameters
//' @return Sample of factor loadings variance, tau
//' @noRd
// [[Rcpp::export]]

arma::vec tau_simulation(List attributes, List hyp_parm, List parm){
  
  int N = attributes["N"];
  int R = attributes["R"];
  int L = attributes["L"];
  mat H = hyp_parm["H"];
  mat B = parm["B"];
  double n_tau = hyp_parm["n.tau"];
  double n_s2_tau = hyp_parm["n.s2.tau"];
  vec tauout(L);
  
  for(int l=0; l < L; l++){
    double n_tau_aux = n_tau + N - l - 1;
    uvec locone = conv_to<uvec>::from(linspace(0,l,l+1));
    uvec loctwo = conv_to<uvec>::from(linspace(l+1,R-1,R-l-1));
    vec bfix(l+1);
    bfix(l) = 1;
    vec b_tau = conv_to<vec>::from(H.submat(loctwo,locone) * inv_sympd(H.submat(locone,locone)) * conv_to<vec>::from(bfix));
    b_tau = conv_to<vec>::from(B(span(l+1,R-1),span(l)) - conv_to<colvec>::from(b_tau));
    mat StateH_tau = H.submat(loctwo,loctwo) - H.submat(loctwo,locone) * inv_sympd(H.submat(locone,locone)) * H.submat(locone,loctwo);
    double n_s2_tau_aux = n_s2_tau + conv_to<double>::from(conv_to<rowvec>::from(b_tau) * StateH_tau * conv_to<colvec>::from(b_tau));
    tauout(l) = 1 / R::rgamma(n_tau_aux/2, 2/n_s2_tau_aux);
    // The line above is not working, gotta fix it!
  }
  
  return tauout;
  
}
