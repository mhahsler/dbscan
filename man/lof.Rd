\name{lof}
\alias{lof}
\alias{LOF}
\title{Local Outlier Factor Score}
\description{
Calculate the Local Outlier Factor (LOF) score for each data point
using a kd-tree to speed up kNN search.
}
\usage{
lof(x, minPts = 5, ...)
}
\arguments{
  \item{x}{ a data matrix or a dist object. }
  \item{minPts}{ number of nearest neighbors used in defining the local neighborhood of a point (includes the point itself). }
  \item{\dots}{ further arguments are passed on to \code{kNN}. Note: \code{lof} uses \code{sort = TRUE} and \code{sort} cannot be specified here. }
}
\details{
LOF compares the local readability density (lrd) of an point to the lrd of its neighbors. A LOF score of approximately 1 indicates that the lrd around the point
is comparable to the lrd of its neighbors and that the point is not an outlier.
Points that have a substantially lower lrd than
their neighbors are considered outliers and produce scores significantly larger than
1.

If a data matrix is specified, then Euclidean distances and fast nearest neighbor search using a kd-tree is used.

\bold{Note on duplicate points:} If there are more than \code{minPts} duplicates of a point in the data, then LOF the local readability distance will be 0 resulting in an undefined LOF score of 0/0. We set LOF in this case to 1 since there is already enough density from the points in the same location to make them not outliers.
The original paper by Breunig et al (2000) assumes that the points are real duplicates and suggests to remove the duplicates before computing LOF. If duplicate points are removed first, then this LOF implementation in \pkg{dbscan} behaves like the one described by Breunig et al.
}
\value{
A numeric vector of length \code{ncol(x)} containing LOF values for
all data points.
}
\references{
Breunig, M., Kriegel, H., Ng, R., and Sander, J. (2000). LOF: identifying
   density-based local outliers. In \emph{ACM Int. Conf. on Management of Data,}
   pages 93-104. \doi{10.1145/335191.335388}
}

\author{
    Michael Hahsler
}
\seealso{
\code{\link{kNN}}, \code{\link{pointdensity}}, \code{\link{glosh}}.
}
\examples{
set.seed(665544)
n <- 100
x <- cbind(
  x=runif(10, 0, 5) + rnorm(n, sd = 0.4),
  y=runif(10, 0, 5) + rnorm(n, sd = 0.4)
  )

### calculate LOF score with a neighborhood of 3 points
lof <- lof(x, minPts = 3)

### distribution of outlier factors
summary(lof)
hist(lof, breaks = 10, main = "LOF (minPts = 3)")

### plot sorted lof. Looks like outliers start arounf a LOF of 2.
plot(sort(lof), type = "l",  main = "LOF (minPts = 3)",
  xlab = "Points sorted by LOF", ylab = "LOF")

### point size is proportional to LOF and mark points with a LOF > 2
plot(x, pch = ".", main = "LOF (minPts = 3)", asp =1)
points(x, cex = (lof-1)*2, pch = 1, col = "red")
text(x[lof>2,], labels = round(lof, 1)[lof>2], pos = 3)
}
\keyword{model}

