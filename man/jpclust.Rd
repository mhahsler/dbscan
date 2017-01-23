\name{jpclust}
\alias{jpclust}
\title{
  Jarvis-Patrick Clustering
}
\description{
  Fast C++ implementation of the Jarvis-Patrick clustering which first
  builds a sparsified neighborhood graph (using the k nearest neighbors) and then places
  two points in the same cluster if at least one is a nearest neighbor of the
  other and they share at least kt nearest neighbors.
}
\usage{
jpclust(x, k, kt, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix/data.frame (Euclidean distance is used) or a precomputed
    dist object.}
  \item{k}{ Neighborhood size.}
  \item{kt}{ threshold for shared nearest neighbors. }
  \item{...}{ additional arguments are passed on to the k
  nearest neighbor search algorithm. See \code{\link{kNN}} for
  details on how to control the search strategy.}
}
\details{
Fast nearest neighbors search with \code{kNN()} is only used if
\code{x} is a matrix. In this case Euclidean distance is used.
}
\value{
  A vector with cluster assignments.
}
\references{
R. A. Jarvis and E. A. Patrick. 1973. Clustering Using a Similarity Measure Based on Shared Near Neighbors. \emph{IEEE Trans. Comput. 22,} 11 (November 1973), 1025-1034.
}

\seealso{
\code{\link{kNN}}
}

\author{
    Michael Hahsler
}

\examples{
data(iris)
iris <- as.matrix(iris[,1:4])

# use a neighborhood of 10 points and require 5 shared neighbors
cl <- jpclust(iris, k = 10, kt = 5)
pairs(iris, col = cl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ clustering }

