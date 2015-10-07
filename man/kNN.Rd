\name{kNN}
\alias{kNN}
\alias{knn}
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
  \item{x}{ a data matrix.}
  \item{k}{ number of neighbors to find. }
    \item{search}{ nearest neighbor search strategy (one of "kdtree"
    or "linear").}
  \item{sort}{ sort the neighbors by distance? }
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
This implementation uses the ANN library (see Mount and Arya, 2010) for
nearest neighbor search. Linear nearest neighbor search can be used. To speed up nearest neighbor search the kd-tree is used.

Note: self-matches are removed!

\code{bucketSize} and \code{splitRule} influence how the kd-tree is built.
\code{approx} uses the approximate nearest neighbor search implemented in ANN.
All nearest neighbors up to a distance of \code{eps}/(1+\code{approx})
will be considered and all with a distance greater than \code{eps} will not
be considered. The other points might be considered. Note that this results in
some actual nearest neighbors being omitted leading to spurious clusters and noise points. However, the algorithm will enjoy a significant speedup. For more details see Mount and Arya (2010).
}
\value{
%%  ~Describe the value returned
A list with the following components:
  \item{dist }{a matrix with distances. }
  \item{id }{a matrix with ids. }
  \item{k }{number of k used. }
%% ...
}
\references{
David M. Mount and Sunil Arya (2010). ANN: A Library for Approximate Nearest Neighbor Searching, \url{https://www.cs.umd.edu/~mount/ANN/}.
}
%\author{
%%  ~~who you are~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{frNN}} for fixed radius nearest neighbors.
}
\examples{
data(iris)
kNN(iris[,-5], k=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
