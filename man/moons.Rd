\name{moons}
\alias{moons}
\docType{data}
\title{
Moons Data
}
\description{
  Contains 100 2-d points, half of which are contained in two moons or "blobs"" (25 points each blob), and the other half in
  asymmetric facing crescent shapes. The three shapes are all linearly separable. 
}
\usage{data("moons")}
\format{
  A data frame with 100 observations on the following 2 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
  }
}
\details{
This data was generated with the following Python commands using the SciKit-Learn library. 
  dontrun{
    import sklearn.datasets as data
    moons, _ = data.make_moons(n_samples=50, noise=0.05)
    blobs, _ = data.make_blobs(n_samples=50, centers=[(-0.75,2.25), (1.0, 2.0)], cluster_std=0.25)
    test_data = np.vstack([moons, blobs])
  }
}
\source{
See the HDBSCAN notebook from github documentation: \url{http://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html}
}
\references{
1. Pedregosa, Fabian, Gael Varoquaux, Alexandre Gramfort, Vincent Michel, Bertrand Thirion, Olivier Grisel, Mathieu Blondel et al. "Scikit-learn: Machine learning in Python." Journal of Machine Learning Research 12, no. Oct (2011): 2825-2830.
}
\examples{
data(moons)
plot(moons, pch=20)
}
\keyword{datasets}
