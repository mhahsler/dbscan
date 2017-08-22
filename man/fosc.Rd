\name{extractFOSC}
\alias{extractFOSC}
\title{
  Framework for Optimal Selection of Clusters
}
\description{
  Generic reimplementation of the Framework for Optimal Selection
  of Clusters (FOSC; Campello et al, 2013). Can be parameterized to perform unsupervised cluster extraction
  through a stability-based measure, or semisupervised cluster extraction through either
  a constraint-based extraction (with a stability-based tiebreaker) or a mixed (weighted)
  constraint and stability-based objective extraction.
}
\usage{
extractFOSC(x, minPts = 2L, constraints = NA, alpha = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a valid hclust object.}
  \item{minPts}{ numeric; Defaults to 2. Only needed if class-less noise is a valid label in the model. }
  \item{constraints}{ Either a list or matrix of pairwise constraints. 
  If not supplied, an unsupervised measure of stability is used for cluster extraction. See details. }
  \item{alpha}{ numeric; weight between [0, 1] for mixed-objective semi-supervised extraction. Defaults to 0. }
}
\details{
Campello et al (2013) suggested a 'Framework for Optimal Selection of Clusters' (FOSC) for 
optimizing arbitrary objective measures that enable local cuts through cluster hierarchies. 
This function implements this framework for hclust objects.

In contrast to global cluster extraction methods (e.g. 'cutting' the tree at a specific height), FOSC traverses a given hierarchy to find a clustering results that is the globally optimal solution to an arbitrary objective measure. As a result of this, multiple subbranches may be 'cut' in the cluster hierarchy to form the final clustering solution. FOSC requires one parameter (minPts) that represents a \emph{minimum} number of points that constitute a valid cluster. Clusters below this threshold are marked as noise, and given the '0' label in the result. Traversing the original cluster tree using minPts creates a new, simplified cluster tree that is then post-processed to maximize the cost function \eqn{J} for each cluster \eqn{C_i}{Ci}, 

\deqn{\max_{\delta_2, \dots, \delta_k} J = \sum\limits_{i=2}^{k} \delta_i S(C_i)}{ J = \sum \delta S(Ci) for all i clusters, }

for some objective measure, \eqn{S}. Currently, the only objective measure supported is the following generic 'stability-based' measure: 
\deqn{ S(C_i) = \sum_{x_j \in C_i}(\frac{1}{h_{min} (x_j, C_i)} - \frac{1}{h_{max} (C_i)}) }{ S(Ci) = \sum (1/Hmin(Xj, Ci) - 1/Hmax(Ci)) for all Xj in Ci.}

\eqn{\delta_i}{\delta} represents an indicator function, which constrains the solution space such that clusters must be
disjoint (cannot assign more than 1 label to each cluster). 

FOSC supports a form of \emph{semisupervised clustering} through \emph{instance-level} constraints (see Wagstaff et al (2002)
for more details). If constraints are given in the call to \code{extractFOSC}, the following alternative objective function is maximized:

\deqn{J = \frac{1}{2n_c}\sum\limits_{j=1}^n \gamma (x_j)}{J = 1/(2 * nc) \sum \gamma(Xj)}

where \eqn{n_c}{nc} is the total number of constraints given and \eqn{\gamma(x_j)}{\gamma(Xj)} represents the number of
constraints involving object \eqn{x_j}{Xj} that are satisfied. In the case of ties (such as solutions where no constraints
were given), the unsupervised solution is used as a tiebreaker. See Campello et al (2013) for more details. FOSC expects the 
pairwise constraints to be passed as either 1) an \eqn{n(n-1)/2} vector of integers representing
the constraints, where 1 represents should-link, -1 represents should-not-link, and 0 represents no preference using the unsupervised solution (see below for examples). Alternatively, if only a few constraints are needed, a named list
representing the (symmetric) adjacency list can be used, where the names correspond to indices of the points in the
original data, and the values correspond to integer vectors of constraints (positive indices for should-link, negative
indices for should-not-link). Again, see the examples section for a demonstration of this. Asymmetric constraints are not supported.

The unsupervised and semisupervised solutions can be weighted with the \eqn{\alpha} parameter to maximize the
following mixed objective function

\deqn{J = \alpha S(C_i) + (1 - \alpha) \gamma(C_i))}{J = \alpha S(Ci) + (1 - \alpha) \gamma(Ci).}

}
\value{
  \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points (if any).}
  \item{hc }{The original hclust object augmented with the n-1 cluster-wide objective scores from the extraction
  encoded in the 'stability', 'constraint', and 'total' named members.}
}
\references{
Campello, Ricardo JGB, et al. "A framework for semi-supervised and unsupervised optimal extraction of clusters from hierarchies." \emph{Data Mining and Knowledge Discovery} 27.3 (2013): 344-371.

Wagstaff, Kiri Lou, and Claire Cardie. \emph{Intelligent clustering with instance-level constraints}. USA: Cornell University, 2002.
}

\seealso{
\code{\link{hdbscan}}, \code{\link[stats]{cutree}}
}

\author{
  Matt Piekenbrock
}

\examples{
data("moons")

## Regular HDBSCAN using stability-based extraction (unsupervised)
cl <- hdbscan(moons, minPts = 5)
cl$cluster

## Constraint-based extraction from the HDBSCAN hierarchy
## (w/ stability-based tiebreaker (semisupervised))
cl_con <- extractFOSC(cl$hc, minPts = 5,
  constraints = list("12" = c(49, -47)))
cl_con$cluster

## Alternative formulation: Constraint-based extraction from the HDBSCAN hierarchy
## (w/ stability-based tiebreaker (semisupervised)) using distance thresholds
dist_moons <- dist(moons)
cl_con2 <- extractFOSC(cl$hc, minPts = 5,
  constraints = ifelse(dist_moons < 0.1, 1L,
                ifelse(dist_moons > 1, -1L, 0L)))

cl_con2$cluster # same as the second example
}
\keyword{ model }
\keyword{ clustering }

