\name{kNNdist}
\alias{kNNdist}
\alias{kNNdistplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate and plot the k-Nearest Neighbor Distance
}
\description{
Fast calculation of the k-nearest neighbor distances in a matrix of points.
The plot can be used to help find a suitable value for the
\code{eps} neighborhood for DBSCAN. Look for the knee in the plot.
}
\usage{
kNNdist(x, k, ...)
kNNdistplot(x, k = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the data set as a matrix or a dist object. }
  \item{k}{ number of nearest neighbors used (use minPoints). }
  \item{...}{further arguments are passed on to \code{kNN}.}
  }
\details{
See \code{\link{kNN}} for a discussion of the kd-tree related parameters.
}
\value{
\code{kNNdist} returns a numeric vector with the distance to its k nearest
neighbor.
}

\author{
    Michael Hahsler
}

\seealso{
\code{\link{kNN}}.
}% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(iris)
iris <- as.matrix(iris[,1:4])

kNNdist(iris, k=4, search="kd")
kNNdistplot(iris, k=4)
## the knee is around a distance of .5

cl <- dbscan(iris, eps = .5, minPts = 4)
pairs(iris, col = cl$cluster+1L)
## Note: black are noise points
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ plot }
