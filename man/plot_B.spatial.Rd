% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_B.spatial.R
\name{plot_B.spatial}
\alias{plot_B.spatial}
\title{Spatial plots of factor loadings}
\usage{
plot_B.spatial(
  Gibbs,
  areapoly,
  burnin = NA,
  permutation = NA,
  main.bool = TRUE,
  layout.dim = NA
)
}
\arguments{
\item{Gibbs}{Result of Gibbs sampler from DIFM function.}

\item{areapoly}{The polygon of the areas. We can obtain this through \code{readOGR} function from \code{sp} package.}

\item{burnin}{Number of burn-in. If not specified, it uses the first tenths as burn-in period.}

\item{permutation}{Permutation of variables. If not specified, no permutation.}

\item{main.bool}{Add title of the plots.}

\item{layout.dim}{Dimension of panel layout for multiple factor loadings. If not specificed, factor loadings plots are layout in one column.}
}
\value{
Factor loadings map plots
}
\description{
The functions builds maps of factor loadings.
}
