\name{optics}
\alias{optics}
\alias{extractDBSCAN}
\alias{extractXi}
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
optics(x, eps = NULL, minPts = 5, ...)

extractDBSCAN(object, eps_cl)
extractXi(object, xi, minimum = FALSE, correctPredecessors = TRUE)

\method{predict}{optics}(object, newdata = NULL, data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix or a dist object.}
  \item{eps}{ upper limit of the size of the epsilon neighborhood. Limiting the
    neighborhood size improves performance and has no or very little impact on
    the ordering as long as it is not set too low. If not specified, the largest minPts-distance in the data set is used which gives the same result as infinity. }
  \item{minPts}{ the parameter is used to identify dense neighborhoods
    and the reachability distance is calculated as the distance to the
    minPts nearest neighbor. Controls the smoothness of the reachability distribution. Default is 5 points.}
  \item{eps_cl}{ Threshold to identify clusters (eps_cl <= eps).}
  \item{xi}{ Steepness threshold to identify clusters hierarchically using
    the Xi method.}
  \item{object}{ an object of class optics. For predict it needs to contain not just an ordering, but also an extracted clustering. }
  \item{data}{ the data set used to create the OPTICS clustering object. }
  \item{newdata}{ new data set for which cluster membership should be predicted. }
  \item{minimum}{ boolean representing whether or not to extract the minimal (non-overlapping)
    clusters in the Xi clustering algorithm. }
  \item{correctPredecessors}{ boolean Correct a common artifacting by pruning the steep up area for
  	points that have predecessors not in the cluster--found by the ELKI framework, see details below. }
  \item{...}{ additional arguments are passed on to fixed-radius
    nearest neighbor search algorithm. See \code{\link{frNN}} for
    details on how to control the search strategy.}
}
\details{
This implementation of OPTICS (Hahsler et al, 2019) implements the original algorithm as described by
Ankerst et al (1999). OPTICS is an ordering algorithm using similar concepts
to DBSCAN. However, for OPTICS
\code{eps} is only an upper limit for the neighborhood size used to reduce
computational complexity. Note that \code{minPts} in OPTICS has a different
effect then in DBSCAN.
It is used to define dense neighborhoods, but since \code{eps} is typically
set rather high, this does not effect the ordering much. However,
it is also used to calculate the reachability distance
and larger values will make the reachability distance plot smoother.

OPTICS linearly orders the data points such that points which are spatially closest become neighbors in the ordering. The closest analog to this ordering is dendrogram
in single-link hierarchical clustering. The algorithm also calculates the
reachability distance for each point. \code{plot()} produces a reachability-plot
which shows each points reachability distance where the points are sorted
by OPTICS. Valleys represent clusters (the deeper the valley, the more dense
  the cluster) and high points indicate points between clusters.

\code{extractDBSCAN} extracts a clustering from an OPTICS ordering that is
similar to what DBSCAN would produce with an eps set to
\code{eps_cl} (see Ankerst et al, 1999).
The only difference to a DBSCAN clustering is that OPTICS is not able to
assign some border points and reports them instead as noise.

\code{extractXi} extract clusters hiearchically
specified in Ankerst et al (1999) based on the steepness of the reachability plot.
One interpretation of the \code{xi}
parameter is that it classifies clusters by change in relative cluster density.
The used algorithm was originally contributed by the ELKI framework
and is explained in Schubert et al (2018),
but contains a set of fixes.

See \code{\link{frNN}} for more information on the parameters related to nearest neighbor search.
}
\value{
  An object of class 'optics' with components:
    \item{eps }{ value of eps parameter. }
    \item{minPts }{ value of minPts parameter. }
    \item{order }{ optics order for the data points in \code{x}. }
    \item{reachdist }{ reachability distance for each data point in \code{x}. }
    \item{coredist }{ core distance for each data point in \code{x}. }

  If \code{extractDBSCAN} was called, then in
      addition the following components are available:
    \item{eps_cl }{ reachability distance for each point in \code{x}. }
    \item{cluster }{ assigned cluster labels in the order of the data points in
      \code{x}. }

  If \code{extractXi} was called, then in
      addition the following components are available:
    \item{xi}{ Steepness threshold\code{x}. }
    \item{cluster }{ assigned cluster labels in the order of the data points in
      \code{x}.}
    \item{clusters_xi }{ data.frame containing the start and end of each cluster found in the OPTICS ordering. }
}
\references{
Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD international conference on Management of data. ACM Press. pp. 49-60. \doi{10.1145/304181.304187}

Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast Density-Based Clustering with R.
  \emph{Journal of Statistical Software}, 91(1), 1-30. \doi{10.18637/jss.v091.i01}

Erich Schubert, Michael Gertz (2018). Improving the Cluster Structure Extracted from OPTICS Plots. Lernen, Wissen, Daten, Analysen (LWDA 2018). pp. 318-329.
}

\author{
    Michael Hahsler and Matthew Piekenbrock
}

\seealso{
\code{\link{frNN}},
\code{\link{as.reachability}}.
}

\examples{
set.seed(2)
n <- 400

x <- cbind(
  x = runif(4, 0, 1) + rnorm(n, sd=0.1),
  y = runif(4, 0, 1) + rnorm(n, sd=0.1)
  )

plot(x, col=rep(1:4, time = 100))

### run OPTICS (Note: we use the default eps calculation)
res <- optics(x, minPts = 10)
res

### get order
res$order

### plot produces a reachability plot
plot(res)

### plot the order of points in the reachability plot
plot(x, col = "grey")
polygon(x[res$order,])

### extract a DBSCAN clustering by cutting the reachability plot at eps_cl
res <- extractDBSCAN(res, eps_cl = .065)
res

plot(res)  ## black is noise
hullplot(x, res)

### re-cut at a higher eps threshold
res <- extractDBSCAN(res, eps_cl = .1)
res
plot(res)
hullplot(x, res)

### extract hierarchical clustering of varying density using the Xi method
res <- extractXi(res, xi = 0.05)
res

plot(res)
hullplot(x, res)

# Xi cluster structure
res$clusters_xi

### use OPTICS on a precomputed distance matrix
d <- dist(x)
res <- optics(d, minPts = 10)
plot(res)
}
\keyword{ model }
\keyword{ clustering }

