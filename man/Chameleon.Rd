\name{Chameleon}
\alias{Chameleon}
\alias{chameleon}
\alias{chameleon_ds4}
\alias{chameleon_ds5}
\alias{chameleon_ds7}
\alias{chameleon_ds8}
\docType{data}
\title{2D Data Sets used for the CHAMELEON Clustering Algorithm}
\description{
Several 2D data sets used to evaluate the CHAMELEON clustering alogrithm
in the paper by Karypis et al (1999).
}
\usage{data(Chameleon)}
\format{
\code{chameleon_ds4}:
  The format is a 8,000 x 2 data.frame.

\code{chameleon_ds5}:
  The format is a 8,000 x 2 data.frame.

\code{chameleon_ds7}:
  The format is a 10,000 x 2 data.frame.

\code{chameleon_ds8}:
  The format is a 8,000 x 2 data.frame.
}
%\details{
%}
\source{
The data was obtained from \url{http://glaros.dtc.umn.edu/gkhome/cluto/cluto/download}
}
\references{
Karypis, G., EH. Han, V. Kumar (1999):
CHAMELEON: A Hierarchical Clustering Algorithm Using Dynamic Modeling,
\emph{IEEE Computer,} \bold{32}(8): 68--75.
}
\examples{
data(Chameleon)

plot(chameleon_ds4, cex=.1)
}
\keyword{datasets}
