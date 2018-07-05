\name{pointdensity}
\alias{pointdensity}
\alias{density}
\title{Calculate Local Density at Each Data Point}
\description{
Calculate the local density at each data point
as either the number of points in the eps-neighborhood (as used in DBSCAN) or
the kernel density estimate (kde) of a uniform kernel.
The function uses a kd-tree for fast fixed-radius nearest neighbor search.
}
\usage{
pointdensity(x, eps, type = "frequency",
  search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix. }
  \item{eps}{ radius of the eps-neighborhood, i.e., bandwidth of
    the uniform kernel).  }
  \item{type}{ \code{"frequency"} or \code{"density"}. should the raw count of points inside the eps-neighborhood
    or the kde be returned. }
  \item{search, bucketSize, splitRule, approx}{ algorithmic parameters for \code{\link{frNN}}. }
}
\details{
DBSCAN estimates the density around a point as the number of points in the
eps-neighborhood of the point (including the query point itself). The kde using a uniform
kernel is just this count divided by \eqn{2 eps n}, where \eqn{n} is the number
of points in \code{x}.

Points with low local density often indicate noise (see e.g., Wishart (1969)
  and Hartigan (1975)).
}
\value{
A vector of the same length as data points (rows) in \code{x} with the
count or density values for each data point.
}
\references{
WISHART, D. (1969), Mode Analysis: A Generalization of Nearest
  Neighbor which Reduces Chaining Effects,
  in Numerical Taxonomy, Ed., A.J. Cole, Academic Press, 282-311.

John A. Hartigan (1975), Clustering Algorithms, John Wiley \& Sons, Inc., New York, NY, USA.
}
\author{Michael Hahsler}
\seealso{
\code{\link{frNN}}.
}
\examples{
set.seed(665544)
n <- 100
x <- cbind(
  x=runif(10, 0, 5) + rnorm(n, sd=0.4),
  y=runif(10, 0, 5) + rnorm(n, sd=0.4)
  )
plot(x)

### calculate density
d <- pointdensity(x, eps = .5, type = "density")

### density distribution
summary(d)
hist(d, breaks = 10)

### point size is proportional to Density
plot(x, pch = 19, main = "Density (eps = .5)", cex = d*5)

### Wishart (1969) single link clustering method
# 1. remove noise with low density
f <- pointdensity(x, eps = .5, type = "frequency")
x_nonoise <- x[f >= 5,]

# 2. use single-linkage on the non-noise points
hc <- hclust(dist(x_nonoise), method = "single")
plot(x, pch = 19, cex = .5)
points(x_nonoise, pch = 19, col= cutree(hc, k = 4)+1L)
}
\keyword{model}

