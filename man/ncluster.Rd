% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncluster.R
\name{ncluster}
\alias{ncluster}
\alias{nnoise}
\alias{nobs}
\title{Number of Clusters, Noise Points, and Observations}
\usage{
ncluster(object, ...)

nnoise(object, ...)
}
\arguments{
\item{object}{a clustering result object containing a \code{cluster} element.}

\item{...}{additional arguments are unused.}
}
\value{
returns the number if clusters or noise points.
}
\description{
Extract the number of clusters or the number of noise points for
a clustering. This function works with any clustering result that
contains a list element named \code{cluster} with a clustering vector. In
addition, \code{nobs} (see \code{\link[stats:nobs]{stats::nobs()}}) is also available to retrieve
the number of clustered points.
}
\examples{
data(iris)
iris <- as.matrix(iris[, 1:4])

res <- dbscan(iris, eps = .7, minPts = 5)
res

ncluster(res)
nnoise(res)
nobs(res)

# the functions also work with kmeans and other clustering algorithms.
cl <- kmeans(iris, centers = 3)
ncluster(cl)
nnoise(cl)
nobs(res)
}
\seealso{
Other clustering functions: 
\code{\link{dbscan}()},
\code{\link{extractFOSC}()},
\code{\link{hdbscan}()},
\code{\link{jpclust}()},
\code{\link{optics}()},
\code{\link{sNNclust}()}
}
\concept{clustering functions}
