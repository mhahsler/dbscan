\name{dbscan}
\alias{dbscan}
\alias{DBSCAN}
%- Also NEED an '\alias' for EACH other topic documented here.
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
  dbscan(x, eps, minPts = 5, weights = NULL,
    borderPoints = TRUE, search = "kdtree", bucketSize = 10,
    splitRule = "suggest", approx = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix or a dist object.}
  \item{eps}{ size of the epsilon neighborhood.}
  \item{minPts}{ number of minimum points in the eps region (for core points).
      Default is 5 points.}
  \item{weights}{ numeric; weights for the data points. Only needed to perform
    weighted clustering. }
  \item{borderPoints}{ logical; should border points be assigned. The default
    is \code{TRUE} for regular DBSCAN. If
    \code{FALSE} then border points are considered noise
    (see DBSCAN* in Campello et al, 2013).}
  \item{search}{ nearest neighbor search strategy (one of "kdtree"
    or "linear", "dist").}
  \item{bucketSize}{ max size of the kd-tree leafs. }
  \item{splitRule}{ rule to split the kd-tree. One of "STD",
    "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST"
    (SL stands for sliding). "SUGGEST" uses ANNs best guess.}
  \item{approx}{ relative error bound for approximate nearest neighbor searching.}
  \item{...}{ additional arguments are currently unused.}
}
\details{
This implementation of DBSCAN implements the original algorithm as described by
Ester et al (1996). DBSCAN estimates the density around each data point by counting the number of points in a user-specified eps-neighborhood and applies a used-specified minPts thresholds to identify core, border and noise points. In a second step, core points are joined into a cluster if they are density-reachable (i.e., there is a chain of core points where one falls inside the eps-neighborhood of the next). Finally, border points are assigned to clusters. The algorithm only needs
parameters \code{eps} and \code{minPts}.

Border points are arbitrarily assigned to clusters in the original algorithm.
DBSCAN* (see Campello et al 2013) treats all border points as noise points. This is implemented with \code{borderPoints = FALSE}.

Setting parameters for DBSCAN:
\code{minPts} is often set to be dimensionality of the data plus one.
The knee in \code{\link{kNNdistplot}} can be used to find suitable values for \code{eps}.

See \code{\link{kNN}} for more information on the other parameters related to
nearest neighbor search.
}
\value{
  A object of class 'dbscan' with the following components:
    \item{eps }{ value of the eps parameter.}
    \item{minPts }{ value of the minPts parameter.}
    \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points.}
  %% ...
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
\code{\link[fpc]{dbscan}} in \pkg{fpc}.
}

\author{
    Michael Hahsler
}

\examples{
data(iris)
iris <- as.matrix(iris[,1:4])

res <- dbscan(iris, eps = .4, minPts = 4)
res

pairs(iris, col = res$cluster + 1L)

## compare with dbscan from package fpc (only if installed)
if (requireNamespace("fpc", quietly = TRUE)) {
  res2 <- fpc::dbscan(iris, eps = .4, MinPts = 4)
  pairs(iris, col = res2$cluster + 1L)

  ## make sure both version produce the same results
  all(res$cluster == res2$cluster)
  }

## find suitable eps parameter (look at knee)
kNNdistplot(iris, k = 4)


## example data from fpc
set.seed(665544)
n <- 100
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

res <- dbscan::dbscan(x, eps = .2, minPts = 4)
res

plot(x, col=res$cluster + 1L)

\dontrun{
## compare speed against fpc version (if microbenchmark is installed)
if (requireNamespace("microbenchmark", quietly = TRUE)) {
  t_dbscan <- microbenchmark::microbenchmark(
    dbscan::dbscan(x, .2, 4), times = 10, unit = "ms")
  t_dbscan_linear <- microbenchmark::microbenchmark(
    dbscan::dbscan(x, .2, 4, search = "linear"), times = 10, unit = "ms")
  t_fpc <- microbenchmark::microbenchmark(
    fpc::dbscan(x, .2, 4), times = 10, unit = "ms")

  rbind(t_fpc, t_dbscan_linear, t_dbscan)

  boxplot(rbind(t_fpc, t_dbscan_linear, t_dbscan),
    names = c("fpc (R)", "dbscan (linear)", "dbscan (kdtree)"),
    main = "Runtime comparison in ms")

  ## speedup of the kd-tree-based version compared to the fpc implementation
  median(t_fpc$time) / median(t_dbscan$time)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ clustering }

