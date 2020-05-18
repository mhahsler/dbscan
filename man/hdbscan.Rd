\name{hdbscan}
\alias{hdbscan}
\alias{HDBSCAN}
\alias{print.hdbscan}
\alias{plot.hdbscan}
\title{
  HDBSCAN
}
\description{
  Fast implementation of the HDBSCAN (Hierarchical DBSCAN) and its related
  algorithms using Rcpp.
}
\usage{
hdbscan(x, minPts, xdist = NULL,
    gen_hdbscan_tree = FALSE,
    gen_simplified_tree = FALSE)

\method{print}{hdbscan}(x, ...)
\method{plot}{hdbscan}(x, scale="suggest",
    gradient=c("yellow", "red"), show_flat = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix or a dist object. }
  \item{minPts}{ integer; Minimum size of clusters. See details. }
  \item{xdist}{ an optional precomputed dist object. May reduce computation time if supplied. }
  \item{gen_hdbscan_tree}{ logical; should the robust single linkage tree be
  explicitly computed. (see cluster tree in Chaudhuri et al, 2010).}
  \item{gen_simplified_tree}{ logical; should the simplified hierarchy be
  explicitly computed. (see Campello et al, 2013). }
  \item{...}{ additional arguments are passed on to the appropriate S3 methods (such as plotting parameters). }
  \item{scale}{ integer; used to scale condensed tree based on the graphics device. Lower scale results in wider trees. }
  \item{gradient}{ character vector; the colors to build the condensed tree coloring with. }
  \item{show_flat}{ logical; whether to draw boxes indicating the most stable clusters. }

}
\details{
This fast implementation of HDBSCAN (Hahsler et al, 2019) computes the hierarchical cluster tree representing density estimates along with the stability-based flat cluster extraction
proposed by Campello et al. (2013). HDBSCAN essentially computes the hierarchy of all DBSCAN* clusterings, and then uses a stability-based extraction method to find optimal cuts in the hierarchy, thus producing a flat solution.

Additional, related algorithms including the "Global-Local Outlier Score from Hierarchies" (GLOSH) (see section 6 of Campello et al., 2015) outlier scores and ability to cluster based on instance-level constraints (see section 5.3 of Campello et al. 2015) are supported. The algorithms only need the parameter \code{minPts}.

Note that \code{minPts} not only acts as a minimum cluster size to detect, but also as a "smoothing" factor of the density estimates implicitly computed from HDBSCAN.
}
\value{
  A object of class 'hdbscan' with the following components:
    \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points.}
    \item{minPts }{ value of the minPts parameter.}
    \item{cluster_scores }{The sum of the stability scores for each salient ('flat') cluster. Corresponds to
    cluster ids given the in 'cluster' member. }
    \item{membership_prob }{The 'probability' or individual stability of a point within its clusters. Between 0 and 1.}
    \item{outlier_scores }{The outlier score (GLOSH) of each point. }
    \item{hc }{An 'hclust' object of the HDBSCAN hierarchy. }
  %% ...
}
\references{
Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast Density-Based Clustering with R.
  \emph{Journal of Statistical Software}, 91(1), 1-30. \doi{10.18637/jss.v091.i01}

Campello RJGB, Moulavi D, Sander J (2013). Density-Based Clustering Based on Hierarchical Density Estimates. \emph{Proceedings of the 17th Pacific-Asia Conference on Knowledge Discovery in Databases, PAKDD 2013,} Lecture Notes in Computer Science 7819, p. 160.

Campello RJGB, Moulavi D, Zimek A, Sander J (2015). Hierarchical density estimates for data clustering, visualization, and outlier detection. \emph{ACM Transactions on Knowledge Discovery from Data (TKDD),} 10(5):1-51. \doi{10.1145/2733381}
}

\seealso{
\code{\link{dbscan}}
}

\author{
    Matt Piekenbrock
}

\examples{
## cluster the moons data set with HDBSCAN
data(moons)

res <- hdbscan(moons, minPts = 5)
res

plot(res)

plot(moons, col = res$cluster + 1L)

## DS3 from Chameleon
data("DS3")

res <- hdbscan(DS3, minPts = 50)
res

## Plot the simplified tree, highlight the most stable clusters
plot(res, show_flat = TRUE)

## Plot the actual clusters
plot(DS3, col=res$cluster+1L, cex = .5)
}
\keyword{ model }
\keyword{ clustering }
\keyword{ hierarchical }

