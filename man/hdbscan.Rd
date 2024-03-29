% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdbscan.R, R/predict.R
\name{hdbscan}
\alias{hdbscan}
\alias{HDBSCAN}
\alias{print.hdbscan}
\alias{plot.hdbscan}
\alias{coredist}
\alias{mrdist}
\alias{predict.hdbscan}
\title{Hierarchical DBSCAN (HDBSCAN)}
\usage{
hdbscan(
  x,
  minPts,
  gen_hdbscan_tree = FALSE,
  gen_simplified_tree = FALSE,
  verbose = FALSE
)

\method{print}{hdbscan}(x, ...)

\method{plot}{hdbscan}(
  x,
  scale = "suggest",
  gradient = c("yellow", "red"),
  show_flat = FALSE,
  ...
)

coredist(x, minPts)

mrdist(x, minPts, coredist = NULL)

\method{predict}{hdbscan}(object, newdata, data, ...)
}
\arguments{
\item{x}{a data matrix (Euclidean distances are used) or a \link{dist} object
calculated with an arbitrary distance metric.}

\item{minPts}{integer; Minimum size of clusters. See details.}

\item{gen_hdbscan_tree}{logical; should the robust single linkage tree be
explicitly computed (see cluster tree in Chaudhuri et al, 2010).}

\item{gen_simplified_tree}{logical; should the simplified hierarchy be
explicitly computed (see Campello et al, 2013).}

\item{verbose}{report progress.}

\item{...}{additional arguments are passed on.}

\item{scale}{integer; used to scale condensed tree based on the graphics
device. Lower scale results in wider trees.}

\item{gradient}{character vector; the colors to build the condensed tree
coloring with.}

\item{show_flat}{logical; whether to draw boxes indicating the most stable
clusters.}

\item{coredist}{numeric vector with precomputed core distances (optional).}

\item{object}{clustering object.}

\item{newdata}{new data points for which the cluster membership should be
predicted.}

\item{data}{the data set used to create the clustering object.}
}
\value{
\code{hdbscan()} returns object of class \code{hdbscan} with the following components:
\item{cluster }{A integer vector with cluster assignments. Zero indicates
noise points.}
\item{minPts }{ value of the \code{minPts} parameter.}
\item{cluster_scores }{The sum of the stability scores for each salient
(flat) cluster. Corresponds to cluster IDs given the in \code{"cluster"} element.
}
\item{membership_prob }{The probability or individual stability of a
point within its clusters. Between 0 and 1.}
\item{outlier_scores }{The GLOSH outlier score of each point. }
\item{hc }{An \link{hclust} object of the HDBSCAN hierarchy. }

\code{coredist()} returns a vector with the core distance for each data point.

\code{mrdist()} returns a \link{dist} object containing pairwise mutual reachability distances.
}
\description{
Fast C++ implementation of the HDBSCAN (Hierarchical DBSCAN) and its related
algorithms.
}
\details{
This fast implementation of HDBSCAN (Campello et al., 2013) computes the
hierarchical cluster tree representing density estimates along with the
stability-based flat cluster extraction. HDBSCAN essentially computes the
hierarchy of all DBSCAN* clusterings, and
then uses a stability-based extraction method to find optimal cuts in the
hierarchy, thus producing a flat solution.

HDBSCAN performs the following steps:
\enumerate{
\item Compute mutual reachability distance mrd between points
(based on distances and core distances).
\item Use mdr as a distance measure to construct a minimum spanning tree.
\item Prune the tree using stability.
\item Extract the clusters.
}

Additional, related algorithms including the "Global-Local Outlier Score
from Hierarchies" (GLOSH; see section 6 of Campello et al., 2015)
is available in function \code{\link[=glosh]{glosh()}}
and the ability to cluster based on instance-level constraints (see
section 5.3 of Campello et al. 2015) are supported. The algorithms only need
the parameter \code{minPts}.

Note that \code{minPts} not only acts as a minimum cluster size to detect,
but also as a "smoothing" factor of the density estimates implicitly
computed from HDBSCAN.

\code{coredist()}: The core distance is defined for each point as
the distance to the \code{MinPts}'s neighbor. It is a density estimate.

\code{mrdist()}: The mutual reachability distance is defined between two points as
\code{mrd(a, b) = max(coredist(a), coredist(b), dist(a, b))}. This distance metric is used by
HDBSCAN. It has the effect of increasing distances in low density areas.

\code{predict()} assigns each new data point to the same cluster as the nearest point
if it is not more than that points core distance away. Otherwise the new point
is classified as a noise point (i.e., cluster ID 0).
}
\examples{
## cluster the moons data set with HDBSCAN
data(moons)

res <- hdbscan(moons, minPts = 5)
res

plot(res)
plot(moons, col = res$cluster + 1L)

## cluster the moons data set with HDBSCAN using Manhattan distances
res <- hdbscan(dist(moons, method = "manhattan"), minPts = 5)
plot(res)
plot(moons, col = res$cluster + 1L)

## DS3 from Chameleon
data("DS3")

res <- hdbscan(DS3, minPts = 50)
res

## Plot the simplified tree, highlight the most stable clusters
plot(res, show_flat = TRUE)

## Plot the actual clusters (noise has cluster id 0 and is shown in black)
plot(DS3, col = res$cluster + 1L, cex = .5)
}
\references{
Campello RJGB, Moulavi D, Sander J (2013). Density-Based Clustering Based on
Hierarchical Density Estimates. Proceedings of the 17th Pacific-Asia
Conference on Knowledge Discovery in Databases, PAKDD 2013, \emph{Lecture Notes
in Computer Science} 7819, p. 160.
\doi{10.1007/978-3-642-37456-2_14}

Campello RJGB, Moulavi D, Zimek A, Sander J (2015). Hierarchical density
estimates for data clustering, visualization, and outlier detection.
\emph{ACM Transactions on Knowledge Discovery from Data (TKDD),} 10(5):1-51.
\doi{10.1145/2733381}
}
\seealso{
Other clustering functions: 
\code{\link{dbscan}()},
\code{\link{extractFOSC}()},
\code{\link{jpclust}()},
\code{\link{optics}()},
\code{\link{sNNclust}()}
}
\author{
Matt Piekenbrock
}
\concept{HDBSCAN functions}
\concept{clustering functions}
\keyword{clustering}
\keyword{hierarchical}
\keyword{model}
