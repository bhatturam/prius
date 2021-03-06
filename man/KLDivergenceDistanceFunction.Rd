% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllFunctions.R
\name{KLDivergenceDistanceFunction}
\alias{KLDivergenceDistanceFunction}
\title{KL Divergence between two vectors.}
\usage{
KLDivergenceDistanceFunction(A, B, extraArgs)
}
\arguments{
\item{A}{A vector of real values}

\item{B}{A vector of real values of the same size as A}

\item{extraArgs}{NA}
}
\value{
A data frame containing the KL Divergence
between A and B as a single value
}
\description{
Computes the KL divergence between two vectors. This is meant
to be used as a distance function for \code{\link{computePathwayScores}}
}

