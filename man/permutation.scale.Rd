% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation.scale.R
\name{permutation.scale}
\alias{permutation.scale}
\title{Permute the dataset by the largest absolute value in each eigenvector, and scale}
\usage{
permutation.scale(data, n.factors, return.scale = FALSE)
}
\arguments{
\item{data}{The dataset}

\item{n.factors}{Number of factors}

\item{return.scale}{Scale data after permutation}
}
\value{
The permuted and standardized dataset, either in matrix or array.
}
\description{
It finds the vector of permutation to permute data by its largest absolute value in each eigenvector. It sets the order by specified number of factors, and the rest is ordered as they were. The data is permuted, and if needed, scaled.
}
