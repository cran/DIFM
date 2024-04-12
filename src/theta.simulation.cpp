#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @title Simulate thetas
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters of the model
//' @param data Dataset of type matrix
//' @param parm Current estimated parameters
//' @return Sample of temporal components, thetas
//' @noRd
// [[Rcpp::export]]

arma::mat theta_simulation(List attributes, List hyp_parm, arma::mat data, List parm){
  
  int N = attributes["N"];
  int R = attributes["R"];
  mat GG = attributes["GG"];
  mat FF = attributes["FF"];
  mat B = parm["B"];
  mat W = parm["W"];
  colvec m0 = hyp_parm["m0"];
  mat C0 = hyp_parm["C0"];
  vec sigma2 = parm["sigma2"];
  int p = GG.n_cols;
  mat FF0 = B * FF;
  mat a(N,p);
  mat m(N,p);
  cube A(p,R,N);
  mat f(N,R);
  mat e(N,R);
  cube RR(p,p,N);
  cube C(p,p,N);
  cube Q(R,R,N);
  
  a.row(0) = conv_to<rowvec>::from(GG * m0);
  RR.slice(0) = GG * C0 * trans(GG) + W;
  f.row(0) = conv_to<rowvec>::from(FF0 * conv_to<colvec>::from(a.row(0)));
  Q.slice(0) = FF0 * RR.slice(0) * trans(FF0) + diagmat(sigma2);
  A.slice(0) = RR.slice(0) * trans(FF0) * inv_sympd(Q.slice(0));
  e.row(0) = conv_to<rowvec>::from(data.row(0) - f.row(0));
  C.slice(0) = inv_sympd(inv_sympd(RR.slice(0)) + trans(FF0) * diagmat(1/sigma2) * FF0);
  m.row(0) = conv_to<rowvec>::from(C.slice(0) * (inv_sympd(RR.slice(0)) * conv_to<colvec>::from(a.row(0)) + trans(FF0) * diagmat(1/sigma2) * conv_to<colvec>::from(data.row(0))));
  
  for(int tt = 1; tt < N; tt++){
    a.row(tt) = conv_to<rowvec>::from(GG * conv_to<colvec>::from(m.row(tt-1)));
    RR.slice(tt) = GG * C.slice(tt-1) * trans(GG) + W;
    f.row(tt) = conv_to<rowvec>::from(FF0 * conv_to<colvec>::from(a.row(tt)));
    Q.slice(tt) = FF0 * RR.slice(tt) * trans(FF0) + diagmat(sigma2);
    A.slice(tt) = RR.slice(tt) * trans(FF0) * inv_sympd(Q.slice(tt));
    e.row(tt) = data.row(tt) - f.row(tt);
    C.slice(tt) = inv_sympd(inv_sympd(RR.slice(tt)) + trans(FF0) * diagmat(1/sigma2) * FF0);
    m.row(tt) = conv_to<rowvec>::from(C.slice(tt) * (inv_sympd(RR.slice(tt)) * conv_to<colvec>::from(a.row(tt)) + trans(FF0) * diagmat(1/sigma2) * conv_to<colvec>::from(data.row(tt))));
  }
  
  mat C_sqrt = sqrtmat_sympd(C.slice(N-1));
  mat thetaout(N,p);
  thetaout.row(N-1) = m.row(N-1) + conv_to<rowvec>::from(C_sqrt * randn(p));
  
  for(int tt = N-2; tt >= 0; tt--){
    mat B = C.slice(tt) * trans(GG) * pinv(RR.slice(tt+1));
    rowvec h = m.row(tt) + conv_to<rowvec>::from(B * conv_to<colvec>::from((thetaout.row(tt+1) - a.row(tt+1))));
    mat H = C.slice(tt) - B * RR.slice(tt+1) * trans(B);
    
    thetaout.row(tt) = h + conv_to<rowvec>::from(sqrtmat_sympd(H) * randn(p));
  }
  
  

//   C.sqrt = matrix.sqrt(C[model.attributes$N,,])
//     theta[model.attributes$N,] <- m[model.attributes$N,] + C.sqrt %*% rnorm(p, mean = 0, sd = 1)
//     
//     for(tt in (model.attributes$N - 1):1)
//     {
//       B <- C[tt,,] %*% t(model.attributes$GG) %*% matrix.Moore.Penrose(R[tt+1,,])
//       h <- m[tt,] + B %*% (theta[tt+1,] - a[tt+1,])
//       H <- C[tt,,] - B %*% R[tt+1,,] %*% t(B)   
//       
//       H.sqrt = matrix.sqrt(H)
//       theta[tt,] <- h + H.sqrt %*% rnorm(p, mean = 0, sd = 1)
//     }
  
  return thetaout;
  
}

