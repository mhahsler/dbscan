\name{sNNclust}
\alias{sNNclust}
\alias{snnclust}
\title{
  Shared Nearest Neighbor Clustering
}
\description{
  Implements the shared nearest neighbor clustering algorithm by Ertoz, Steinbach
  and Kumar.
}
\usage{
sNNclust(x, k, minstr, minlinks, ...)
}
\arguments{
  \item{x}{a data matrix/data.frame (Euclidean distance is used), a precomputed
    dist object or a kNN object created with \code{kNN()}.}
  \item{k}{ Neighborhood size for nearest neighbor sparsification to create the
    shared NN graph. Link strength is the number of shared NN. }
  \item{minstr}{ link strength threshold for the shared NN graph. }
  \item{minlinks}{ threshold for representative points. I.e., the number
    of links in the shared NN graph.}
  \item{...}{ additional arguments are passed on to the k
  nearest neighbor search algorithm. See \code{\link{kNN}} for
  details on how to control the search strategy.}
}
%\details{
%}
\value{
  A vector with cluster assignments.
}
\references{
Levent Ertoz, Michael Steinbach, Vipin Kumar, A New Shared Nearest Neighbor Clustering Algorithm and its Applications, \emph{Workshop on Clustering High Dimensional Data and its Applications at 2nd SIAM International Conference on Data Mining,} 2002.
}

\seealso{
\code{\link{jpclust}}
}

\author{
    Michael Hahsler
}

\examples{
data(iris)
iris <- as.matrix(iris[,1:4]) ### only use numeric attributes

# Out of 20 NN 10 have to be shared (minstr) to create a link in the sNN graph.
# A point needs a least 5 strong links to be a representative (not noise).
# Noise points have cluster id 0 and are shown in black.
cl <- sNNclust(iris, k = 20, minstr = 10, minlinks = 5)
pairs(iris, col = cl + 1L)
}
\keyword{ model }
\keyword{ clustering }

