\name{sNN}
\alias{sNN}
\alias{snn}
\alias{sort.sNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shared Nearest Neighbors}
\description{
Calculates the number of shared nearest neighbors.
}
\usage{
sNN(x, k, kt = NULL, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix, a dist object or a kNN object. }
  \item{k}{ number of neighbors to consider to calculate the shared nearest neighbors. }
  \item{kt}{ threshold on the number of shared nearest neighbors graph. Edges
    are only preserved if \code{kt} or more neighbors are shared. }
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
plot(sNN(nn, kt = 3), x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
