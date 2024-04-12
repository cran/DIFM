
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

std::vector<double> dinvgamma(arma::vec x, double alpha, double beta, bool logtrans){
  
  int xn = x.n_elem;
  arma::vec Yvec(xn);
  NumericVector Yvec_nv(xn);
  for(int i = 0; i < xn; i++){
    Yvec(i) = pow(beta, alpha) / tgamma(alpha) * pow(1./x(i), alpha + 1) * exp(-beta / x(i));
  }
  
  if(logtrans){Yvec = log(Yvec);}
  
  std::vector<double> Yvec_std(Yvec.begin(), Yvec.end());
  
  for(int i = 0; i < xn; i++){
    Yvec_nv[i] = Yvec(i);
  }
  
  return Yvec_std;
  
}


double dinvwishart(arma::mat Sigma, double nu, arma::mat S, bool logtrans) {
  
  int k = Sigma.n_rows;
  double dens = 0, gamsum = 0;
  for (int i = 0; i < k; i++) {
    gamsum += lgamma((nu + 1.0 - i - 1.0) / 2.0);
  }
  dens = -((nu * k) / 2.0) * log(2.0) - ((k * (k - 1.0)) / 4.0) * 
    log(M_PI) - gamsum + (nu / 2.0) * log(arma::det(S)) - 
    ((nu + k + 1.0) / 2.0) * log(arma::det(Sigma)) - 
    0.5 * arma::trace(S * arma::inv(Sigma));
  if (logtrans == false) dens = exp(dens);
  return dens;
}
// Copied from LaplacesDemonCPP 



// [[Rcpp::depends(RcppArmadillo)]];

//' @title  Marginal predictive density 
//' @description  It calculates the marginal density (Lewis and Raftery, 1997) from the DIFM sample using C++.
//' @param  data The dataset
//' @param  attributes Model attributes generated from \code{difm.model.attributes}.
//' @param  hyp_parm Hyperparameters generated from \code{difm.hyp.parm}.
//' @param  Gibbs Result of Gibbs sampler from DIFM function.
//' @param  burnin Burn-in period. If not specified, one tenths of the iterations will be the burn-in period.
//' @param  verbose  Print out the process.
//'
//' @return  A list of 4 items: Laplace-Metropolis predictive density of the given DIFM, integrated likelihood, the maximum of the predictive densities and determinant of the covariance matrix of the parameters.
//' 
//' @export
// [[Rcpp::export]]

List marginal_d_cpp(arma::mat data, List attributes, List hyp_parm, List Gibbs, int burnin = -1, bool verbose = true){
  
  double N = attributes["N"];
  double R = attributes["R"];
  double L = attributes["L"];
  mat GG = attributes["GG"];
  mat FF = attributes["FF"];
  colvec m0 = hyp_parm["m0"];
  mat C0 = hyp_parm["C0"];
  double p = GG.n_cols;
  double n_tau = hyp_parm["n.tau"];
  double n_s2_tau = hyp_parm["n.s2.tau"];
  double n_w = hyp_parm["n.w"];
  double n_sigma = hyp_parm["n.sigma"];
  double n_s2_sigma = hyp_parm["n.s2.sigma"];
  mat Psi = hyp_parm["Psi"];
  mat H = hyp_parm["H"];
  mat Hplus = hyp_parm["Hplus"];

  cube B = Gibbs["B"];
  cube theta = Gibbs["theta"];
  mat sigma2 = Gibbs["sigma2"];
  mat tau = Gibbs["tau"];
  cube W = Gibbs["W"];
  double n_size = sigma2.n_rows;

  if(burnin == -1){
    burnin = floor(n_size/10);
  }
  mat tot_result_mat(n_size - burnin, N);
  mat intgr_like_mat(n_size - burnin, N);
  
  double npars = R*L - L*(L+1)/2 + R + 2*pow(L, 2) + L + L;
  mat parmat(n_size - burnin, npars);
  int parmat_i = 0;
  for(int i1 = 0; i1 < L; i1++){
    for(int i2 = i1 + 1; i2 < R; i2++){
      vec Badd = B(span(burnin, n_size - 1), span(i2), span(i1));
      parmat(span(0, n_size - burnin - 1), span(parmat_i)) = conv_to<colvec>::from(Badd);
      parmat_i += 1;
    }
  }
  for(int i1 = 0; i1 < p; i1++){
    for(int i2 = i1; i2 < p; i2++){
      vec Wadd = W(span(burnin, n_size - 1), span(i1), span(i2));
      parmat(span(0, n_size - burnin - 1), span(parmat_i)) = conv_to<colvec>::from(Wadd);
      parmat_i += 1;
    }
  }
  for(int i = 0; i < R; i++){
    vec sigma2add = sigma2(span(burnin, n_size - 1), span(i));
    parmat(span(0, n_size - burnin - 1), span(parmat_i)) = conv_to<colvec>::from(sigma2add);
    parmat_i += 1;
  }
  for(int i = 0; i < L; i++){
    vec tauadd = tau(span(burnin, n_size - 1), span(i));
    parmat(span(0, n_size - burnin - 1), span(parmat_i)) = conv_to<colvec>::from(tauadd);
    parmat_i += 1;
  }

  mat parcov = cov(parmat);
  vec par_eigval;
  mat par_eigvec;
  eig_sym(par_eigval, par_eigvec, parcov);
  double detparcov = accu(log(par_eigval));

  for(int g = burnin; g < n_size; g++){
    vec tot_result(N, arma::fill::value(npars/2*log(2*M_PI)/N + 0.5*detparcov/N));
    mat W0 = W(span(g), span(0, p-1), span(0, p-1));
    mat B0 = B(span(g), span(0, R-1), span(0, L-1));
    if(L == 1){B0 = trans(B0);}
    vec sigma20 = conv_to<vec>::from(sigma2.row(g));
    vec tau0 = conv_to<vec>::from(tau.row(g));
    mat V0 = diagmat(sigma2.row(g));
    mat FF0 = B0 * FF;
    mat a(N,p);
    mat m(N,p);
    cube A(p,R,N);
    mat f(N,R);
    mat e(N,R);
    cube RR(p,p,N);
    cube C(p,p,N);
    cube Q(R,R,N);

    a.row(0) = conv_to<rowvec>::from(GG * m0);
    RR.slice(0) = GG * C0 * trans(GG) + W0;
    f.row(0) = conv_to<rowvec>::from(FF0 * conv_to<colvec>::from(a.row(0)));
    Q.slice(0) = FF0 * RR.slice(0) * trans(FF0) + V0;
    A.slice(0) = RR.slice(0) * trans(FF0) * inv_sympd(Q.slice(0));
    e.row(0) = conv_to<rowvec>::from(data.row(0) - f.row(0));
    C.slice(0) = inv_sympd(inv_sympd(RR.slice(0)) + trans(FF0) * inv_sympd(V0) * FF0);
    m.row(0) = conv_to<rowvec>::from(C.slice(0) * (inv_sympd(RR.slice(0)) * conv_to<colvec>::from(a.row(0)) + trans(FF0) * inv_sympd(V0) * conv_to<colvec>::from(data.row(0))));

    double log_det_Q0 = arma::log_det_sympd(Q.slice(0)) + conv_to<double>::from(e.row(0)*inv_sympd(Q.slice(0))*conv_to<colvec>::from(e.row(0)));
    tot_result(0) = tot_result(0) -(R/2)*log(2*M_PI) - 0.5*log_det_Q0;
    
    for(int tt = 1; tt < N; tt++){
      a.row(tt) = conv_to<rowvec>::from(GG * conv_to<colvec>::from(m.row(tt-1)));
      RR.slice(tt) = GG * C.slice(tt-1) * trans(GG) + W0;
      f.row(tt) = conv_to<rowvec>::from(FF0 * conv_to<colvec>::from(a.row(tt)));
      Q.slice(tt) = FF0 * RR.slice(tt) * trans(FF0) + V0;
      A.slice(tt) = RR.slice(tt) * trans(FF0) * inv_sympd(Q.slice(tt));
      e.row(tt) = data.row(tt) - f.row(tt);
      C.slice(tt) = inv_sympd(inv_sympd(RR.slice(tt)) + trans(FF0) * inv_sympd(V0) * FF0);
      m.row(tt) = conv_to<rowvec>::from(C.slice(tt) * (inv_sympd(RR.slice(tt)) * conv_to<colvec>::from(a.row(tt)) + trans(FF0) * inv_sympd(V0) * conv_to<colvec>::from(data.row(tt))));

      double log_det_Qtt = log_det_sympd(Q.slice(tt)) + conv_to<double>::from(e.row(tt)*inv_sympd(Q.slice(tt))*conv_to<colvec>::from(e.row(tt)));
      tot_result(tt) = tot_result(tt) -(R/2)*log(2*M_PI) - 0.5*log_det_Qtt;
    }

    double mar1 = accu(as<arma::vec>(wrap(dinvgamma(sigma20, n_sigma/2, n_s2_sigma/2, true))));
    double mar2 = dinvwishart(W0, n_w + N - 1, Psi, true);
    double mar3 = accu(as<arma::vec>(wrap(dinvgamma(tau0, n_tau/2, n_s2_tau/2, true))));
    double mar4 = 0;

    for(int i = 0; i < L; i++){
      uvec locone = linspace<uvec>(0, i, i+1);
      uvec loctwo = linspace<uvec>(i+1, R-1, R-i-1);
      colvec mufix = zeros<colvec>(i+1);
      mufix(i) = 1;
      colvec b_b = zeros<colvec>(R-i-1);
      b_b = b_b -Hplus.submat(loctwo, locone) * inv_sympd(Hplus(locone, locone)) * mufix;
      mat b_Hplus = Hplus(loctwo, loctwo) - Hplus(loctwo, locone) * inv_sympd(Hplus(locone, locone)) * Hplus(locone, loctwo);
      vec b_Hplus_eigval;
      mat b_Hplus_eigvec;
      eig_sym(b_Hplus_eigval, b_Hplus_eigvec, b_Hplus);
      vec Bnow = conv_to<vec>::from(B0.col(i));
      double log_b_mar = -(R-2-i)/2*log(2*M_PI) - 0.5*accu(log(b_Hplus_eigval(span(1, R-i-2)))) - 0.5* conv_to<double>::from(trans(conv_to<colvec>::from(Bnow(span(i+1, R-1))) - b_b) * pinv(b_Hplus) * (conv_to<colvec>::from(Bnow(span(i+1, R-1))) - b_b));
      mar4 += log_b_mar;
    }

    vec mar1rep(N, arma::fill::value(mar1/N));
    vec mar2rep(N, arma::fill::value(mar2/N));
    vec mar3rep(N, arma::fill::value(mar3/N));
    vec mar4rep(N, arma::fill::value(mar4/N));
    tot_result_mat.row(g-burnin) = conv_to<rowvec>::from(tot_result + mar1rep + mar2rep + mar3rep + mar4rep);
    intgr_like_mat.row(g-burnin) = conv_to<rowvec>::from(tot_result);

    if(verbose && (g+1) % 100 == 0){ Rcout << "Current step: " << g+1 << "th" << std::endl;}
  }

  vec tot_result_rsum = conv_to<vec>::from(arma::sum(tot_result_mat, 1));
  double Maximum = max(tot_result_rsum);
  List result_list = List::create(Named("Total") = tot_result_mat, Named("IL") = intgr_like_mat, Named("Maximum") = Maximum, Named("detparcov") = detparcov);
  
  return result_list;
  
}
