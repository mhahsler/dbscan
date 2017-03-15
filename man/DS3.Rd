\name{DS3}
\alias{DS3}
\docType{data}
\title{
DS3: Spatial data with arbitrary shapes
}
\description{
  Contains 8000 2-d points, with 6 "natural" looking shapes, all of which have an sinesoid-like shape that intersects 
  with each cluster. 
}
\usage{data("DS3")}
\format{
  A data frame with 8000 observations on the following 2 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
  }
}
\details{
Originally used as a benchmark data set for the Chameleon clustering algorithm[1] to illustrate the a data set
containing arbitrarily shaped spatial data surrounded by both noise and artifacts. 
}
\source{
Obtained at https://cs.joensuu.fi/sipu/datasets/
}
\references{
1. Karypis, George, Eui-Hong Han, and Vipin Kumar. "Chameleon: Hierarchical clustering using dynamic modeling." Computer 32, no. 8 (1999): 68-75.
}
\examples{
data(DS3)
plot(DS3, pch=20, cex=0.25)
}
\keyword{datasets}
