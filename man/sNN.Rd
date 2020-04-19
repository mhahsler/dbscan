\name{sNN}
\alias{sNN}
\alias{snn}
\alias{sort.sNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shared Nearest Neighbors}
\description{
Calculates the number of shared nearest neighbors, the shared nearest neighbor similarity and
creates a shared nearest neighbors graph.
}
\usage{
sNN(x, k, kt = NULL, jp = FALSE, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix, a dist object or a kNN object. }
  \item{k}{ number of neighbors to consider to calculate the shared nearest neighbors. }
  \item{kt}{ minimum threshold on the number of shared nearest neighbors to build the shared
    nearest neighbor graph. Edges
    are only preserved if \code{kt} or more neighbors are shared. }
  \item{jp}{ use the definition by Javis and Patrick (1973), where shared neighbors are only counted between points that are in each other's neighborhood, otherwise 0 is returned. If FALSE, then the number of shared neighbors is returned, even if the points are not neighbors. }
  \item{search}{ nearest neighbor search strategy (one of "kdtree", "linear"
      or "dist").}
  \item{sort}{ sort by the number of shared nearest neighbors? Note that this is expensive and
      \code{sort = FALSE} is much faster. sNN objects can be sorted using
      \code{sort()}.}
  \item{bucketSize}{ max size of the kd-tree leafs. }
  \item{splitRule}{ rule to split the kd-tree. One of "STD",
      "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST"
      (SL stands for sliding). "SUGGEST" uses ANNs best guess.}
  \item{approx}{ use approximate nearest neighbors. All NN up to a distance of
    a factor of 1+\code{approx} eps may be used. Some actual NN may be
    omitted leading to spurious clusters and noise points.
    However, the algorithm will enjoy a significant speedup. }
}
\details{
The number of shared nearest neighbors is the intersection of the
kNN neighborhood of two points. Note: that each point is considered to be part
of its own kNN neighborhood. The range for the shared nearest neighbors is
[0,k].

Javis and Patrick (1973) use the shared nearest neighbor graph for clustering. They only count shared neighbors between points that are in each other's kNN neighborhood.
}
\value{
An object of class sNN (subclass of kNN and NN) containing a list with the following components:
  \item{id }{a matrix with ids. }
  \item{dist }{a matrix with the distances. }
  \item{shared }{a matrix with the number of shared nearest neighbors. }
  \item{k }{number of k used. }
%% ...
}
\seealso{
\code{\link{NN}} and \code{\link{kNN}} for k nearest neighbors.
}

\author{
    Michael Hahsler
}

\references{
R. A. Jarvis and E. A. Patrick. 1973. Clustering Using a Similarity Measure Based on Shared Near Neighbors. \emph{IEEE Trans. Comput. 22,} 11 (November 1973), 1025-1034.
}

\examples{
data(iris)
x <- iris[, -5]

# finding kNN and add the number of shared nearest neighbors.
k <- 5
nn <- sNN(x, k = k)
nn

# shared nearest neighbor distribution
table(as.vector(nn$shared))

# explore neighborhood of point 10
i <- 10
nn$shared[i,]

plot(nn, x)

# apply a threshold to create a sNN graph with edges
# if more than 3 neighbors are shared.
nn_3 <- sNN(nn, kt = 3)
plot(nn_3, x)

# get an adjacency list for the shared nearest neighbor graph
adjacencylist(nn_3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
