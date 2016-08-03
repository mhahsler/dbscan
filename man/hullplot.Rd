\name{hullplot}
\alias{hullplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Convex Hulls of Clusters
}
\description{
    This function produces a two-dimensional scatter plot with added
    convex hulls for clusters.
}
\usage{
hullplot(x, cl, col = NULL, cex = 0.5, hull_lwd = 1, hull_lty = 1,
	 main = "Convex Cluster Hulls", solid=FALSE, alpha = .2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data matrix (only the first two columns are used to produce
      the scatter plot).}
  \item{cl}{a clustering. Either a numeric cluster assignment vector or a clustering object (a list with an element named \code{cluster}).}
  \item{col}{ colors used for clusters. Defaults to the standard palette.
    The first color is used for noise/unassigned points (cluster id 0). }
  \item{cex}{ expansion factor for symbols. }
  \item{hull_lwd, hull_lty}{ line width and line type used for the convex hull.}
  \item{main}{ main title.}
  \item{solid, alpha}{ draw filled polygons instead of lines for the convex hulls.
    alpha controls the alpha shading. }
  \item{...}{ additional arguments passed on to plot. }
}
%\details{
%}
%\value{}

\author{
    Michael Hahsler
}

%\seealso{
%}

\examples{
set.seed(2)
n <- 400

x <- cbind(
  x = runif(4, 0, 1) + rnorm(n, sd=0.1),
  y = runif(4, 0, 1) + rnorm(n, sd=0.1)
  )
cl <- rep(1:4, time = 100)

### original data with true clustering
### (also with differnt symbols, just the hulls or polygons)
hullplot(x, cl, main = "True clusters")
hullplot(x, cl, main = "True clusters", pch = cl)
hullplot(x, cl, main = "True clusters", pch = NA)
hullplot(x, cl, main = "True clusters", solid = TRUE, hull_lwd = 0)


### run some clustering algorithms and plot the resutls
db <- dbscan(x, eps = .07, minPts = 10)
hullplot(x, db, main = "DBSCAN")

op <- optics(x, eps = 10, minPts = 10, eps_cl = 0.07)
hullplot(x, op, main = "OPTICS")

opXi <- opticsXi(op, xi = 0.05)
hullplot(x, opXi, main = "OPTICSXi")

km <- kmeans(x, centers = 4)
hullplot(x, km, main = "k-means")

hc <- cutree(hclust(dist(x)), k = 4)
hullplot(x, hc, main = "Hierarchical")
}
\keyword{ plot }
\keyword{ clustering }

