% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllFunctions.R
\name{computePersonalizationVectors}
\alias{computePersonalizationVectors}
\title{Compute the personalization vectors based on experimental data based on a
score function}
\usage{
computePersonalizationVectors(experimentalData, PPIGraph, defaults,
  scoreFunction, scoreFunctionExtraArgs)
}
\arguments{
\item{experimentalData}{experimentalData  A data frame with as many rows as
HGNC Gene Symbols, with the columns defined by the output of testFunction.
The output of \code{\link{runTestOnData}}.}

\item{PPIGraph}{PPIGraph A list containing igraph_object - the igraph object
and vertex_map - a frame containing node attributes currently indegree and
outdegree. The output of \code{\link{createIGraphObject}}}

\item{defaults}{A data frame corresspoding to the default experimental values
for genes missing from the experiment. This must be of the same format as
the output of the scoreFunction.}

\item{scoreFunction}{A score function that accepts the experimentalData and
the PPI Graph to return the normal and disease personalization vectors. Foe
example see \code{\link{outdegreeNormalizedFCScoreFunction}} or
\code{\link{outdegreeNormalizedFCOneMinusPScoreFunction}}}

\item{scoreFunctionExtraArgs}{A list of extra arguments that are to be passed
to the score function if any.}
}
\value{
A data frame with two columns npv - normal personalization vector and
  dpv - disease personalization vector. There are as many rows as number of
  HGNC Symbols in PPIGraph
}
\description{
This function computes the personalization vectors from the experiment data
based on an input score function based on the test parameters.
}
\examples{
\dontrun{
gds = getGEO('GDS3837',AnnotGPL = TRUE, getGPL = TRUE,destdir = "/tmp")
gpl = getGEO(Meta(gds)$platform,destdir = "/tmp")
eset = GDS2eSet(gds, GPL=gpl,do.log2=FALSE)
numPairs = dim(pData(eset))[1]/2
expressionDataGEO=importExpressionDataGEO(eset,gpl,defaultProbeSelector,probeCombinerMean)
experimentDataGEO=runTestOnData(expressionDataGEO,1:60,61:120,logPairedTTestFunction)
personalizationVectorsGEO=computePersonalizationVectors(experimentDataGEO,reactomePPIGraph,data.frame(foldChange=1,pValue=1),outdegreeNormalizedFCOneMinusPScoreFunction)
}
}

