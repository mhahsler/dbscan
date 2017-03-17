\name{jpclust}
\alias{jpclust}
\title{
  Jarvis-Patrick Clustering
}
\description{
  Fast C++ implementation of the Jarvis-Patrick clustering which first
  builds a shared nearest neighbor graph (k nearest neighbor sparsification)
  and then places
  two points in the same cluster if they are in each other's nearest neighbor
  list and they share at least kt nearest neighbors.
}
\usage{
jpclust(x, k, kt, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix/data.frame (Euclidean distance is used), a precomputed
    dist object or a kNN object created with \code{kNN()}.}
  \item{k}{ Neighborhood size for nearest neighbor sparsification. If \code{x} is
a kNN object then \code{k} may be missing. }
  \item{kt}{ threshold on the number of shared nearest neighbors (including the
    points themselves) to form clusters. }
  \item{...}{ additional arguments are passed on to the k
  nearest neighbor search algorithm. See \code{\link{kNN}} for
  details on how to control the search strategy.}
}
\details{
Note: Following the original paper, the shared nearest neighbor list is
constructed as the k neighbors plus the point itself (as neighbor zero).
Therefore, the threshold \code{kt} can be in the range [1, k].

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
data("Chameleon")

# use a shared neighborhood of 20 points and require 12 shared neighbors
cl <- jpclust(chameleon_ds4, k = 20, kt = 12)
plot(chameleon_ds4, col = cl, cex = .5)
# Note how the sine wave points chain clusters together.

# use a precomputed kNN object instead of the original data.
nn <- kNN(chameleon_ds4, k = 30)
nn

cl <- jpclust(nn, k = 20, kt = 12)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ clustering }

