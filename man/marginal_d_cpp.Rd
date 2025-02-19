% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{marginal_d_cpp}
\alias{marginal_d_cpp}
\title{Marginal predictive density}
\usage{
marginal_d_cpp(data, attributes, hyp_parm, Gibbs, burnin = -1L, verbose = TRUE)
}
\arguments{
\item{data}{The dataset}

\item{attributes}{Model attributes generated from \code{difm.model.attributes}.}

\item{hyp_parm}{Hyperparameters generated from \code{difm.hyp.parm}.}

\item{Gibbs}{Result of Gibbs sampler from DIFM function.}

\item{burnin}{Burn-in period. If not specified, one tenths of the iterations will be the burn-in period.}

\item{verbose}{Print out the process.}
}
\value{
A list of 4 items: Laplace-Metropolis predictive density of the given DIFM, integrated likelihood, the maximum of the predictive densities and determinant of the covariance matrix of the parameters.
}
\description{
It calculates the marginal density (Lewis and Raftery, 1997) from the DIFM sample using C++.
}
