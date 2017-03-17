\name{NN}
\alias{adjacencylist}
\alias{adjacencylist.kNN}
\alias{adjacencylist.frNN}
\alias{plot.NN}
\alias{NN}
\alias{nn}
\title{Nearest Neighbors Auxiliary Functions}
\description{
Helper functions for nearest neighbors (\link{kNN} and \link{frNN}).
}
\usage{
adjacencylist(x, ...)
\S3method{plot}{NN}(x, data, main = NULL, ...)
}
\arguments{
  \item{x}{ a nearest neighbor object (of class kNN or frNN). }
  \item{...}{ further arguments are currently ignored. }
  \item{data}{ data with the coordinates for plotting. }
  \item{main}{ main title for the plot. }
}
%\details{
%}
\value{
%%  ~Describe the value returned
\code{adjacencylist} returns
a list with one element for each original data point. Each element contains
the row ids of the nearest neighbors. The adjacency list can be used to create
for example a graph object.
}
\seealso{
\code{\link{frNN}} and \code{\link{kNN}}.
}

\author{
    Michael Hahsler
}

\examples{
data(iris)
x <- iris[, -5]

# finding kNN directly in data (using a kd-tree)
nn <- kNN(x, k=5)
nn

# plot the kNN where NN are shown as line conecting points.
plot(nn, x)

# show the first few elements of the adjacency list
head(adjacencylist(nn))

# create a graph and find connected components
library("igraph")
g <- graph_from_adj_list(adjacencylist(nn))
comp <- components(g)
plot(x, col = comp$membership)

# detect clusters (communities) with the label propagation algorithm
cl <- membership(cluster_label_prop(g))
plot(x, col = cl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
