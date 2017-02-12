\name{sNNclust}
\alias{sNNclust}
\alias{snnclust}
\title{Shared Nearest Neighbor Clustering}
\description{
  Implements the shared nearest neighbor clustering algorithm by Ertoz, Steinbach
  and Kumar.
}
\usage{
sNNclust(x, k, eps, minPts, borderPoints = TRUE, ...)
}
\arguments{
  \item{x}{a data matrix/data.frame (Euclidean distance is used), a precomputed
    dist object or a kNN object created with \code{kNN()}.}
  \item{k}{ Neighborhood size for nearest neighbor sparsification to create the
    shared NN graph. }
  \item{eps}{ Two objects are only reachable from each other if they share
    at least \code{eps} nearest neighbors. }
  \item{minPts}{ minimum number of points that share at least \code{eps}
    nearest neighbors for a point to be considered a core points.  }
  \item{borderPoints}{ should borderPoints be assigned to clusters like in
    DBSCAN?}
  \item{...}{ additional arguments are passed on to the k
    nearest neighbor search algorithm. See \code{\link{kNN}} for
    details on how to control the search strategy.}
}
\details{
Algorithm:

1) Constructs a
shared nearest neighbor graph for a given k. The edge weights are the number of shared k nearest neighbors (in the range of [0, k]).

2) Find each points SNN density, i.e., the number of points which have a similarity of eps or greater.

3) Find the core points, i.e., all points that have an SNN density
greater than MinPts.

4) Form clusters from the core points and assign border points (i.e., non-core points which share at least eps neighbors with a core point).

Note that steps 2-4 are DBSCAN and that \code{eps} is used on a similarity (the number of shared neighbors) and not on a distance like in DBSCAN.
}
\value{
  A vector with cluster assignments.
}
\references{
Levent Ertoz, Michael Steinbach, Vipin Kumar, Finding Clusters of Different Sizes, Shapes, and Densities in Noisy, High Dimensional Data, \emph{SIAM International Conference on Data Mining,} 2003, 47-59.
}

\seealso{
\code{\link{jpclust}}
}

\author{
    Michael Hahsler
}

\examples{
data("Chameleon")

# Out of k = 20 NN 7 (eps) have to be shared to create a link in the sNN graph.
# A point needs a least 16 (minPts) links in the sNN graph to be a core point.
# Noise points have cluster id 0 and are shown in black.
cl <- sNNclust(chameleon_ds4, k = 20, eps = 7, minPts = 16)
plot(chameleon_ds4, col = cl + 1L, cex = .5)
}
\keyword{ model }
\keyword{ clustering }

