\name{kNN}
\alias{kNN}
\alias{knn}
\alias{sort.kNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the k Nearest Neighbors}
\description{
This function uses a kd-tree to find all k nearest neighbors
in a data matrix (including distances) fast.
}
\usage{
kNN(x, k, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix, a dist object or a kNN object. }
  \item{k}{ number of neighbors to find. }
    \item{search}{ nearest neighbor search strategy (one of "kdtree", "linear"
      or "dist").}
  \item{sort}{ sort the neighbors by distance? Note that this is expensive and
      \code{sort = FALSE} is much faster. kNN objects can be sorted using
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
\code{search} controls if a kd-tree or linear search (both implemented in
the ANN library; see Mount and Arya, 2010). Note, that these implementations cannot handle NAs. \code{search="dist"} precomputes Euclidean distances first using R.
NAs are handled, but the resulting distance matrix cannot contain NAs. To use other distance measures, a precomputed distance matrix can be
provided as \code{x} (\code{search} is ignored).

\code{bucketSize} and \code{splitRule} influence how the kd-tree is built.
\code{approx} uses the approximate nearest neighbor search implemented in ANN.
All nearest neighbors up to a distance of \code{eps}/(1+\code{approx})
will be considered and all with a distance greater than \code{eps} will not
be considered. The other points might be considered. Note that this results in
some actual nearest neighbors being omitted leading to spurious clusters and noise points. However, the algorithm will enjoy a significant speedup. For more details see Mount and Arya (2010).

\emph{Note:} self-matches are removed!
}
\value{
%%  ~Describe the value returned
An object of class kNN containing a list with the following components:
  \item{dist }{a matrix with distances. }
  \item{id }{a matrix with ids. }
  \item{k }{number of k used. }
%% ...
}
\references{
David M. Mount and Sunil Arya (2010). ANN: A Library for Approximate Nearest Neighbor Searching, \url{http://www.cs.umd.edu/~mount/ANN/}.
}
\seealso{
\code{\link{NN}} and \code{\link{frNN}} for fixed radius nearest neighbors.
}

\author{
    Michael Hahsler
}

\examples{
data(iris)
x <- iris[, -5]

# finding kNN directly in data (using a kd-tree)
nn <- kNN(x, k=5)
nn

# explore neighborhood of point 10
i <- 10
nn$id[i,]
plot(x, col = ifelse(1:nrow(iris) \%in\% nn$id[i,], "red", "black"))

# Visualize the 5 nearest neighbors
plot(nn, x)

# Visualize a reduced 2-NN graph
plot(kNN(nn, k = 2), x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
