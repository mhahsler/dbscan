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
LOF compares the local density of an point to the local densities
of its neighbors. Points that have a substantially lower density than
their neighbors are considered outliers.
A LOF score of approximately 1 indicates that density around the point
is comparable to its neighbors. Scores significantly larger than
1 indicate outliers.

Note: If there are more than \code{k} duplicate points in the data, then lof
can become \code{NaN} caused by an infinite local density.
In this case we set lof to 1. The paper by Breunig et al (2000) suggests a different method of removing all duplicate points first.
}
\value{
A numeric vector of length \code{ncol(x)} containing LOF values for
all data points.
}
\references{
Breunig, M., Kriegel, H., Ng, R., and Sander, J. (2000). LOF: identifying
   density-based local outliers. In \emph{ACM Int. Conf. on Management of Data,}
   pages 93-104.
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

