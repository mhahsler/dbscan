\name{reachability}
\alias{reachability}
\alias{as.reachability}
\alias{as.reachability.optics}
\alias{as.reachability.dendrogram}
\alias{as.dendrogram}
\alias{as.dendrogram.reachability}
\alias{plot.reachability}
\alias{print.reachability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Density Reachability Structures
}
\description{
  Class "reachability" provides general functions for representing various hierarchical representations as "reachability plots", as originally defined by Ankerst et al (1999). Methods include fast implementations of the conversion algorithms introduced by Sanders et al (2003) to convert between dendrograms and reachability plot objects.
}
\usage{
as.reachability(object, ...)

\method{as.reachability}{optics}(object, ...)

\method{as.reachability}{dendrogram}(object, ...)

\method{as.dendrogram}{reachability}(object, ...)

\method{print}{reachability}(x, ...)

\method{plot}{reachability}(x, order_labels = FALSE,
  xlab = "Order", ylab = "Reachability dist.",
  main = "Reachability Plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{any R object that can be made into one of class "reachability", such as an object of class "optics" or "dendrogram".}
  \item{x}{object of class "reachability".}
  \item{order_labels}{ whether to plot text labels for each points reachability distance.}
  \item{xlab}{ x-axis label, defaults to "Order".}
  \item{ylab}{ y-axis label, defaults to "Reachability dist.".}
  \item{main}{ Title of the plot, defaults to "Reachability Plot".}
  \item{...}{ graphical parameters, or arguments for other methods. }
}
\details{
Dendrograms are a popular visualization tool for representing hierarchical relationships. In agglomerative clustering, dendrograms can be constructed using a variety of linkage criterion (such as single or complete linkage), many of which are frequently used to 1) visualize the density-based relationships in the data or 2) extract cluster labels from the data the dendrogram represents.

The original ordering algorithm OPTICS as described by Ankerst et al (1999) introduced the notion of
2-dimensional representation of so-called "density-reachability" that was shown to be useful for data
visualization. This representation was shown to essentially convey the same information as the more
traditional dendrogram structure by Sanders et al (2003).

Different hierarchical representations, such as dendrograms or reachability plots, may be preferable depending on the context. In smaller datasets, cluster memberships may be more easily identifiable through a dendrogram representation, particularly is the user is already familiar with tree-like representations. For larger datasets however, a reachability plot may be prefered for visualizing macro-level density relationships.

The central idea behind a reachability plot is that the ordering in which points are plotted identifies
underlying hierarchical density representation. OPTICS linearly orders the data points such that points which are spatially closest become neighbors in the ordering. Valleys represent clusters, which can be represented hierarchically. Although the ordering is crucial to the structure of the reachability plot, its important to note that OPTICS, like DBSCAN, is not entirely deterministic and, just like the dendrogram, isomorphisms may exist

A variety of cluster extraction methods have been proposed using reachability plots. Because both cluster
extraction depend directly on the ordering OPTICS produces, they are part of the optics interface. Noneless,
reachability plots can be created directly from other types of linkage trees, and vice versa.

See \code{\link{optics}} for more information on how OPTICS is formulated. \code{\link{extractDBSCAN}} and \code{\link{extractXi}} are the two cluster extraction methods presented in the original OPTICS publication.
}
\value{
  An object of class 'reachability' with components:
    \item{order }{ order to use for the data points in \code{x}. }
    \item{reachdist }{ reachability distance for each data point in \code{x}. }
}
\references{
Ankerst, M., M. M. Breunig, H.-P. Kriegel, J. Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. \emph{ACM SIGMOD international conference on Management of data.} ACM Press. pp. 49--60.

Sander, J., X. Qin, Z. Lu, N. Niu, and A. Kovarsky (2003). Automatic extraction of clusters from hierarchical clustering representations. \emph{Pacific-Asia Conference on Knowledge Discovery and Data Mining.} Springer Berlin Heidelberg.
}

\author{
    Matthew Piekenbrock
}
\seealso{
  \code{\link[stats]{dendrogram}},
  \code{\link{optics}},
  \code{\link{extractDBSCAN}},
  \code{\link{extractXi}},
  \code{\link{hclust}}.
}

\examples{
set.seed(2)
n <- 20

x <- cbind(
  x = runif(4, 0, 1) + rnorm(n, sd=0.1),
  y = runif(4, 0, 1) + rnorm(n, sd=0.1)
)

plot(x, xlim=range(x), ylim=c(min(x)-sd(x), max(x)+sd(x)), pch=20)
text(x = x, labels = 1:nrow(x), pos=3)

### run OPTICS
res <- optics(x, eps = 10,  minPts = 2)
res

### plot produces a reachability plot
plot(res)

### Extract reachability components from OPTICS
reach <- as.reachability(res)
reach

### plot still produces a reachability plot; points ids
### (rows in the original data) can be displayed with order_labels = TRUE
plot(reach, order_labels = TRUE)

### Reachability objects can be directly converted to dendrograms
dend <- as.dendrogram(reach)
dend
plot(dend)

### A dendrogram can be converted back into a reachability object
plot(as.reachability(dend))


}
\keyword{ model }
\keyword{ clustering }
\keyword{ hierarchical clustering }

