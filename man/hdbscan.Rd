\name{hdbscan}
\alias{hdbscan}
\alias{HDBSCAN}
\alias{plot.hdbscan}
\title{
  HDBSCAN
}
\description{
  Fast implementation of the HDBSCAN (Hierarchical DBSCAN) and its related 
  algorithms using Rcpp. 
}
\usage{
hdbscan(x, minPts, xdist = NULL, gen_rsl_tree = FALSE, gen_condensed_tree = FALSE, ...)

\method{plot}{hdbscan}(object, scale="suggest", gradient=c("yellow", "red"), show_flat = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix or matrix-coercible object. }
  \item{minPts}{ integer; Minimum size of clusters. See details. }
  \item{xdist}{ dist object (euclidean required); if available. }
  \item{gen_rsl_tree}{ logical; should the robust single linkage tree be 
  explicitly computed. (see cluster tree in Chaudhuri et al, 2010).}
  \item{gen_condensed_tree}{ logical; should the simplified hierarchy be 
  explicitly computed. (see Campello et al, 2013). }
  \item{...}{ additional arguments are passed on to the appropriate S3 methods when specialized, 
  otherwise unused. }
  \item{object}{ a HDBSCAN clustering object. }
  \item{scale}{ integer; used to scale condensed tree based on the graphics device. Lower scale results in wider trees. }
  \item{gradient}{ character vector; the colors to build the condensed tree coloring with. }
  \item{show_flat}{ logical; whether to draw boxes indicating the most stable clusters. }  

}
\details{
Computes the hierarchical "cluster tree" representing density estimates along with the stability-based 'flat' cluster extraction 
proposed by Campello et al. (2013). HDBSCAN essentially computes the hierarchy of all DBSCAN* (again, see Campello et al. 2013) clusterings, and then uses a stability-based extraction method to find optimal 'cuts' in the hierarchy, thus producing a flat solution. Additional, related algorithms including the "Global-Local Outlier Score from Hierarchies" (GLOSH) (see section 6 of Campello et al. 2015) outlier scores and ability to cluster based on instance-level constraints (see section 5.3 of Campello et al. 2015) are supported. The algorithm only need the parameter \code{minPts}.

Note that \code{minPts} not only acts as a minimum cluster size to detect, but also as a "smoothing" factor of the density 
estimates implicitly computed from HDBSCAN. 
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
Martin Ester, Hans-Peter Kriegel, Joerg Sander, Xiaowei Xu (1996). A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise. Institute for Computer Science, University of Munich. \emph{Proceedings of 2nd International Conference on Knowledge Discovery and Data Mining (KDD-96).}

Campello, R. J. G. B.; Moulavi, D.; Sander, J. (2013). Density-Based Clustering
Based on Hierarchical Density Estimates. \emph{Proceedings of the 17th
Pacific-Asia Conference on Knowledge Discovery
in Databases, PAKDD 2013,} Lecture Notes in Computer Science 7819, p. 160.

Campello, Ricardo JGB, et al. "Hierarchical density estimates for data clustering, visualization, and outlier detection." ACM Transactions on Knowledge Discovery from Data (TKDD) 10.1 (2015): 5.

Chaudhuri, Kamalika, and Sanjoy Dasgupta. "Rates of convergence for the cluster tree." Advances in Neural Information Processing Systems. 2010.
}

\seealso{
\code{\link{dbscan}}
}

\author{
    Michael Hahsler, 
    Matt Piekenbrock
}

\examples{
data(moons)

res <- hdbscan(moons, minPts = 5)
res

plot(moons, col = res$cluster + 1L)

## example data from fpc
set.seed(665544)
n <- 100
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

res <- hdbscan(x, minPts = 3)
res

## Plot the simplified tree, highlight the most stable clusters 
plot(res, show_flat = T)

## Plot the actual clusters
plot(x, col=res$cluster+1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ clustering }

