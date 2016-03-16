\name{optics}
\alias{optics}
\alias{optics_cut}
\alias{predict.optics}
\alias{OPTICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  OPTICS
}
\description{
  Implementation of the OPTICS
  (Ordering points to identify the clustering structure) clustering
  algorithm using a kd-tree.
}
\usage{
  optics(x, eps, minPts = 5, eps_cl,
    search = "kdtree", bucketSize = 10,
    splitRule = "suggest", approx = 0)

  optics_cut(x, eps_cl)

  \method{predict}{optics}(object, data, newdata = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix or a distance matrix.}
  \item{eps}{ upper limit of the size of the epsilon neighborhood.}
  \item{minPts}{ number of minimum points in the eps region (for core points).
      Default is 5 points.}
  \item{eps_cl}{ Threshold to identify clusters (eps_cl <= eps).}
  \item{search}{ nearest neighbor search strategy (one of "kdtree",
    "linear" search, or precomputed "dist") Using precomputed distances
    is better for high-dimensional, but smaller data sets.}
  \item{bucketSize}{ max size of the kd-tree leafs. }
  \item{splitRule}{ rule to split the kd-tree. One of "STD",
      "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST"
      (SL stands for sliding). "SUGGEST" uses ANNs best guess.}
  \item{approx}{ relative error bound for approximate nearest neighbor searching.}
  \item{object}{ an OPTICS object that contains a clustering. }
  \item{data}{ the data set used to create the OPTICS clustering object. }
  \item{newdata}{ new data set for which cluster membership should be predicted. }
  \item{...}{ additional arguments are currently unused.}
}
\details{
This implementation of OPTICS implements the original algorithm as described by
Ankers et al (1999). OPTICS is similar to DBSCAN, however, for OPTICS
\code{eps} is only an upper limit for the neighborhood size used to reduce
computational complexity. Similar to DBSCAN, \code{minPts} is often set to be dimensionality of the data plus one.

OPTICS linearly orders the data points such that points which are spatially closest become neighbors in the ordering. The closest analog to this ordering is dendrogram
in single-link hierarchical clustering. The algorithm also calculates the
reachability distance for each point. \code{plot()} produces a reachability-plot
which shows each points reachability distance where the points are sorted
by OPTICS. Valleys represent clusters (the deeper the valley, the more dense
  the cluster) and high points indicate points between clusters.

If \code{eps_cl} is specified, then an algorithm to extract clusters
(see Ankers et al, 1999) is used. That is, it internally calls \code{optics_cut}
to extract the clustering.
The resulting clustering is similar to what DBSCAN would
produce. The only difference is that OPTICS is not able to assign some border
points and reports them instead as noise.

See \code{\link{kNN}} for more information on the other parameters related to
nearest neighbor search.
}
\value{
  An object of class 'optics' with components:
    \item{eps }{ value of eps parameter. }
    \item{minPts }{ value of minPts parameter. }
    \item{order }{ optics order for the data points in \code{x}. }
    \item{reachdist }{ reachability distance for each data point in \code{x}. }
    \item{coredist }{ core distance for each data point in \code{x}. }

  If \code{eps_cl} was specified or \code{optics_cut} was called, then in
      addition the following components are available:
    \item{eps_cl }{ reachability distance for each point in \code{x}. }
    \item{cluster }{ assigned cluster labels in the order of the data points in
      \code{x}. }
}
\references{
Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD international conference on Management of data. ACM Press. pp. 49--60.
}

\author{
    Michael Hahsler
}

\seealso{
\code{\link[fpc]{dbscan}} in \pkg{fpc}.
}

\examples{
set.seed(2)
n <- 400

x <- cbind(
  x = runif(4, 0, 1) + rnorm(n, sd=0.1),
  y = runif(4, 0, 1) + rnorm(n, sd=0.1)
  )

plot(x, col=rep(1:4, time = 100))

### run OPTICS
res <- optics(x, eps = 1,  minPts = 10)
res

### get order
res$order

### plot produces a reachability plot
plot(res)

### identify clusters by cutting the reachability plot (black is noise)
res <- optics_cut(res, eps_cl = .065)
res

plot(res)
plot(x, col = res$cluster+1L)

### re-cutting at a higher eps threshold
res <- optics_cut(res, eps_cl = .1)
res
plot(res)
plot(x, col = res$cluster+1L)

### using OPTICS on a precomputed distance matrix
d <- dist(x)
res <- optics(x, eps = 1, minPts = 10)
plot(res)
}
\keyword{ model }
\keyword{ clustering }

