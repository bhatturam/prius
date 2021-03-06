% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllFunctions.R
\name{logPairedTTestFunction}
\alias{logPairedTTestFunction}
\title{Paired t-test based test function}
\usage{
logPairedTTestFunction(nData, dData, hgncSymbol, extraArgs)
}
\arguments{
\item{nData}{A vector with numSamples/2  normal expression values for a gene}

\item{dData}{A vector with numSamples/2  disease expression values for a gene}

\item{hgncSymbol}{The HGNC Symbol of the gene}

\item{extraArgs}{#' @param extraArgs A list of named extra arguments, the
only name supported is exponent (assumed to be 1 if not provided for non
log transformed data). An exponent not equal to 1 assumes that the
expression data has been log transformed with base exponent}
}
\value{
A single row data frame with two columns foldChange and pValue
}
\description{
Perform a paired t-test on expression data for a single gene
It is meant to be used with \code{\link{runTestOnData}}
This function is used to obtain a p-value on the fold change between normal
and disease samples for a single gene. If the data is already logarithm
transformed in this function, then the exponent must be passed in the list
extraArgs
}

