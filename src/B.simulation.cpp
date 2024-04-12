
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title Simulate Bs
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters of the model
//' @param data Dataset of type matrix
//' @param parm Current estimated parameters
//' @return Sample of factor loadings, B
//' @noRd
// [[Rcpp::export]]
 
arma::mat B_simulation(List attributes, List hyp_parm, arma::mat data, List parm){
  
  int N = attributes["N"];
  int L = attributes["L"];
  int R = attributes["R"];
  arma::mat H = hyp_parm["H"];
  arma::mat Hplus = hyp_parm["Hplus"];
  arma::mat X = parm["X"];
  arma::vec tau = parm["tau"];
  arma::vec sigma2 = parm["sigma2"];
  arma::mat Vinv = diagmat(1/sigma2);
  arma::mat Bvar(L*R, L*R);
  arma::mat Bvarinv(L*R, L*R);
  arma::colvec mufix(accu(linspace(1,L,L)));
  arma::vec kseq(accu(linspace(1,L,L)));
  arma::colvec kseqnot = linspace(0,R*L-1,R*L);
  arma::colvec mu_b(R*L);
  arma::mat xstar;
  int kfill = 0;
  arma::vec mufixaux = cumsum(linspace(1,L,L));
  for(int i = 0; i < L; i++){
    mufix(mufixaux(i) - 1) = 1;
    Bvar(span(R*i,R*(i+1)-1),span(R*i,R*(i+1)-1)) = Hplus * tau(i);
    Bvarinv(span(R*i,R*(i+1)-1),span(R*i,R*(i+1)-1)) = H / tau(i);
    for(int i2 = 0; i2 <= i; i2++){
      kseq(kfill) = R*(i) + i2; // In terms of R code, it should be R*(i) + i2 + 1
      kfill += 1;
    }
  }
  uvec ksequ = conv_to<uvec>::from(kseq);
  kseqnot.shed_rows(ksequ);
  uvec kseqnotu = conv_to<uvec>::from(kseqnot);
  mat Sigbinv = Bvarinv;
  for(int i = 0; i < N; i++){
    xstar = kron(X.row(i),eye(R,R));
    Sigbinv += trans(xstar) * Vinv * xstar;
    mu_b += trans(xstar) * Vinv * conv_to<colvec>::from(data.row(i));
  }
  mat Sigb = inv_sympd(Sigbinv);
  mu_b = Sigb * mu_b;
  mat mu_b_star = mu_b.rows(kseqnotu) + Sigb.submat(kseqnotu,ksequ) * inv_sympd(Sigb.submat(ksequ,ksequ)) * (mufix - mu_b.rows(ksequ));
  mat Sig_b_star = Sigb.submat(kseqnotu, kseqnotu) - Sigb.submat(kseqnotu, ksequ) * pinv(Sigb.submat(ksequ, ksequ)) * Sigb.submat(ksequ, kseqnotu);
  vec b_aux = mu_b_star + sqrtmat_sympd(Sig_b_star) * randn(R*L - accu(linspace(1,L,L)));
  mat Bout(R,L);
  Bout.diag().ones();
  for(int i = 0; i < L; i++){
    Bout(span(i+1,R-1),span(i)) = b_aux(span(R*i - accu(linspace(0,i,i+1)), R*(i+1) - accu(linspace(1,i+1,i+1)) - 1));
  }
  
  
  
  return Bout;
  
}
 
