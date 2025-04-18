% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getTopMarkers}
\alias{getTopMarkers}
\title{Extract Top Markers from Seurat FindAllMarkers Results}
\usage{
getTopMarkers(mks, n = 5)
}
\arguments{
\item{mks}{Data frame output from Seurat FindAllMarkers, including 'gene' and 'cluster' columns.}

\item{n}{Number of genes to extract per cluster. Default is 5.}
}
\value{
A vector of unique top markers from each cluster.
}
\description{
Extract Top Markers from Seurat FindAllMarkers Results
}
