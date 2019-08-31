\name{dbscan}
\alias{dbscan}
\alias{DBSCAN}
\alias{predict.dbscan_fast}
\title{
  DBSCAN
}
\description{
  Fast reimplementation of the DBSCAN (Density-based spatial
  clustering of applications with noise) clustering algorithm
  using a kd-tree. The implementation is
  significantly faster and can work with larger data sets then dbscan in
  \pkg{fpc}.
}
\usage{
dbscan(x, eps, minPts = 5, weights = NULL, borderPoints = TRUE, ...)

\method{predict}{dbscan_fast}(object, newdata = NULL, data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix or a dist object. Alternatively, a frNN object
    with fixed-radius nearest neighbors can also be specified (see Example section). In this case
    eps can be missing and will be taken form the frNN object.}
  \item{eps}{ size of the epsilon neighborhood.}
  \item{minPts}{ number of minimum points in the eps region (for core points).
      Default is 5 points.}
  \item{weights}{ numeric; weights for the data points. Only needed to perform
    weighted clustering. }
  \item{borderPoints}{ logical; should border points be assigned. The default
    is \code{TRUE} for regular DBSCAN. If
    \code{FALSE} then border points are considered noise
    (see DBSCAN* in Campello et al, 2013).}
  \item{object}{ a DBSCAN clustering object. }
  \item{data}{ the data set used to create the DBSCAN clustering object. }
  \item{newdata}{ new data set for which cluster membership should be predicted. }
  \item{...}{ additional arguments are passed on to fixed-radius
  nearest neighbor search algorithm. See \code{\link{frNN}} for
  details on how to control the search strategy.}
}
\details{
\emph{Note:} use \code{dbscan::dbscan} to call this
implementation when you also use package \pkg{fpc}.

This implementation of DBSCAN implements the original algorithm as described by
Ester et al (1996). DBSCAN estimates the density around each data point by counting the number of points in a user-specified eps-neighborhood and applies a used-specified minPts thresholds to identify core, border and noise points. In a second step, core points are joined into a cluster if they are density-reachable (i.e., there is a chain of core points where one falls inside the eps-neighborhood of the next). Finally, border points are assigned to clusters. The algorithm only needs
parameters \code{eps} and \code{minPts}.

Border points are arbitrarily assigned to clusters in the original algorithm.
DBSCAN* (see Campello et al 2013) treats all border points as noise points. This is implemented with \code{borderPoints = FALSE}.

Setting parameters for DBSCAN:
\code{minPts} is often set to be dimensionality of the data plus one or higher.
The knee in \code{\link{kNNdistplot}} can be used to find suitable values for \code{eps}.

See \code{\link{frNN}} for more information on the parameters related to nearest neighbor search.

A precomputed frNN object can be supplied as \code{x}. In this case
\code{eps} does not need to be specified. This option us useful
for large data sets, where a sparse distance matrix is available.
See \code{\link{frNN}} how to create frNN objects.

\code{predict} can be used to predict cluster memberships for new data points.
A point is considered a member of a cluster if it is within the eps
neighborhood of a member of the cluster. Points which cannot be assigned to a
cluster will be reported as members of the noise cluster 0.
}
\value{
  An object of class 'dbscan_fast' with the following components:
    \item{eps }{ value of the eps parameter.}
    \item{minPts }{ value of the minPts parameter.}
    \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points.}
}
\references{
Martin Ester, Hans-Peter Kriegel, Joerg Sander, Xiaowei Xu (1996). A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise. Institute for Computer Science, University of Munich. \emph{Proceedings of 2nd International Conference on Knowledge Discovery and Data Mining (KDD-96).}

Campello, R. J. G. B.; Moulavi, D.; Sander, J. (2013). Density-Based Clustering
Based on Hierarchical Density Estimates. \emph{Proceedings of the 17th
Pacific-Asia Conference on Knowledge Discovery
in Databases, PAKDD 2013,} Lecture Notes in Computer Science 7819, p. 160.
}

\seealso{
\code{\link{kNNdistplot}},
\code{\link{frNN}},
\code{\link[fpc]{dbscan}} in \pkg{fpc}.
}

\author{
    Michael Hahsler
}

\examples{
## Example 1: use dbscan on the iris data set
data(iris)
iris <- as.matrix(iris[,1:4])

## find suitable eps parameter using a k-NN plot for k = dim + 1
## Look for the knee!
kNNdistplot(iris, k = 5)
abline(h=.5, col = "red", lty=2)

res <- dbscan(iris, eps = .5, minPts = 5)
res

pairs(iris, col = res$cluster + 1L)

## use precomputed frNN
fr <- frNN(iris, eps = .5)
dbscan(fr, minPts = 5)

## Example 2: use data from fpc
set.seed(665544)
n <- 100
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

res <- dbscan(x, eps = .3, minPts = 3)
res

## plot clusters and add noise (cluster 0) as crosses.
plot(x, col=res$cluster)
points(x[res$cluster==0,], pch = 3, col = "grey")

hullplot(x, res)

## predict cluster membership for new data points
## (Note: 0 means it is predicted as noise)
newdata <- x[1:5,] + rnorm(10, 0, .2)
predict(res, newdata, data = x)

## compare speed against fpc version (if microbenchmark is installed)
## Note: we use dbscan::dbscan to make sure that we do now run the
## implementation in fpc.
\dontrun{
if (requireNamespace("fpc", quietly = TRUE) &&
    requireNamespace("microbenchmark", quietly = TRUE)) {
  t_dbscan <- microbenchmark::microbenchmark(
    dbscan::dbscan(x, .3, 3), times = 10, unit = "ms")
  t_dbscan_linear <- microbenchmark::microbenchmark(
    dbscan::dbscan(x, .3, 3, search = "linear"), times = 10, unit = "ms")
  t_dbscan_dist <- microbenchmark::microbenchmark(
    dbscan::dbscan(x, .3, 3, search = "dist"), times = 10, unit = "ms")
  t_fpc <- microbenchmark::microbenchmark(
    fpc::dbscan(x, .3, 3), times = 10, unit = "ms")

  r <- rbind(t_fpc, t_dbscan_dist, t_dbscan_linear, t_dbscan)
  r

  boxplot(r,
    names = c('fpc', 'dbscan (dist)', 'dbscan (linear)', 'dbscan (kdtree)'),
    main = "Runtime comparison in ms")

  ## speedup of the kd-tree-based version compared to the fpc implementation
  median(t_fpc$time) / median(t_dbscan$time)
}}

## Example 3: manually create a frNN object for dbscan (dbscan only needs ids and eps)
nn <- structure(list(ids = list(c(2,3), c(1,3), c(1,2,3), c(3,5), c(4,5)), eps = 1),
  class =  c("NN", "frNN"))
nn
dbscan(nn, minPts = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ clustering }

