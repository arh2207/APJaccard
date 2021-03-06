% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/APJaccard.R
\name{JaccardSimTransform}
\alias{JaccardSimTransform}
\title{A function that calculates a transformed Jaccard similarity matrix for the cells in a Seurat object.
tanh with k=280 performed the best (outperformed Louvain) in one example
with an annotated data set.}
\usage{
JaccardSimTransform(data.dist, k = 5, transformation = "default")
}
\arguments{
\item{data.dist}{a distance matrix on which to use FindNeighbors}

\item{k}{an integer representing the number of neighbors which will be found for each cell}

\item{transformation}{a character vector indicating the type of transformation the similarity matrix will go through before being returned. By default the similarity matrix will be returned with no transformation. Other possible values include log (returns matrix which has been log transformed), 'inverse' returns a shifted inverse of the matrix, 'stretch' returns a linear scaling of the matrix, 'tan' returns the tangent of matrix values, and 'tanh' returns the hyperbolic tangent of matrix values}
}
\value{
an NxN matrix where N is the number of cells, and the matrix contains the transformed Jaccard values for all pairs of cells
}
\description{
A function that calculates a transformed Jaccard similarity matrix for the cells in a Seurat object.
tanh with k=280 performed the best (outperformed Louvain) in one example
with an annotated data set.
}
