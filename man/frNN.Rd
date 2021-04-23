\name{frNN}
\alias{frNN}
\alias{frnn}
\alias{sort.frNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the Fixed Radius Nearest Neighbors}
\description{
This function uses a kd-tree to find the fixed radius nearest neighbors
(including distances) fast.
}
\usage{
frNN(x, eps, query = NULL, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix, a dist object or a frNN object.}
  \item{eps}{ neighbors radius. }
  \item{query}{ a data matrix with the points to query. If query is not specified, the NN for all the points in \code{x} is returned. If query is specified then \code{x} needs to be a data matrix. }
  \item{sort}{ sort the neighbors by distance? This is expensive and can be done later using \code{sort()}. }
  \item{search}{ nearest neighbor search strategy (one of "kdtree"
    or "linear", "dist").}
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
\emph{Note:} self-matches are not returned!

If \code{x} is specified as a data matrix, then Euclidean distances an
fast nearest neighbor lookup using a kd-tree are used.
See \code{\link{kNN}} for details on the parameters for the kd-tree.

To create a frNN object from scratch, you need to supply at least the
elements \code{id} with a list of integer vectors with
the nearest neighbor ids
for each point and \code{eps} (see below).
}
\value{
%%  ~Describe the value returned
An object of class frNN (subclass of NN) containing a list with the following components:
  \item{id }{a list of integer vectors. Each vector contains the ids
  of the fixed radius nearest neighbors. }
  \item{dist }{a list with distances (same structure as \code{ids}). }
  \item{eps }{ eps used. }
%% ...
}
\seealso{
\code{\link{NN}} and \code{\link{kNN}} for k nearest neighbor search.
}
\references{
David M. Mount and Sunil Arya (2010). ANN: A Library for Approximate Nearest Neighbor Searching, \url{http://www.cs.umd.edu/~mount/ANN/}.
}

\author{
    Michael Hahsler
}

\examples{
data(iris)
x <- iris[, -5]

# Example 1: Find fixed radius nearest neighbors for each point
nn <- frNN(x, eps=.5)

# Number of neighbors
hist(sapply(adjacencylist(nn), length),
  xlab = "k", main="Number of Neighbors",
  sub = paste("Neighborhood size eps =", nn$eps))

# Explore neighbors of point i = 10
i <- 10
nn$id[[i]]
nn$dist[[i]]
plot(x, col = ifelse(1:nrow(iris) \%in\% nn$id[[i]], "red", "black"))

# get an adjacency list
head(adjacencylist(nn))

# plot the fixed radius neighbors (and then reduced to a radius of .3)
plot(nn, x)
plot(frNN(nn, eps = .3), x)

## Example 2: find fixed-radius NN for query points
q <- x[c(1,100),]
nn <- frNN(x, eps = .5, query = q)

plot(nn, x, col = "grey")
points(q, pch = 3, lwd = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
