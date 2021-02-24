\name{lof}
\alias{lof}
\alias{LOF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Local Outlier Factor Score}
\description{
Calculate the Local Outlier Factor (LOF) score for each data point
using a kd-tree to speed up kNN search.
}
\usage{
lof(x, k = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix or a dist object. }
  \item{k}{ size of the neighborhood. }
  \item{\dots}{ further arguments are passed on to \code{kNN}. }
}
\details{
LOF compares the local reachability density (lrd) of an point to the lrd of its neighbors. A LOF score of approximately 1 indicates that the lrd around the point
is comparable to the lrd of its neighbors and that the point is not an outlier.
Points that have a substantially lower lrd than
their neighbors are considered outliers and produce scores significantly larger than
1.

\bold{Note on duplicate points:} If there are more than \code{k} duplicates of a point in the data, then LOF
can become \code{NaN} caused by an lrd of 0 restulting leading to LOF of0/0.
For some applications, the duplicate points should be in the data (e.g., they are the result of rounding the values for several points).
We set LOF to 1 in this case since there is already enought density from the points in the same location to make them not outliers.
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
  x=runif(10, 0, 5) + rnorm(n, sd=0.4),
  y=runif(10, 0, 5) + rnorm(n, sd=0.4)
  )

### calculate LOF score
lof <- lof(x, k=3)

### distribution of outlier factors
summary(lof)
hist(lof, breaks=10)

### point size is proportional to LOF
plot(x, pch = ".", main = "LOF (k=3)")
points(x, cex = (lof-1)*3, pch = 1, col="red")
text(x[lof>2,], labels = round(lof, 1)[lof>2], pos = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}

