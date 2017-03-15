\name{extractFOSC}
\alias{extractFOSC}
\title{
  Framework for Optimal Selection of Clusters
}
\description{
  Generic reimplementation of the Framework for Optimal Selection 
  of Clusters (FOSC). Can be parameterized to perform unsupervised cluster extraction 
  through a stability-based measure, or semisupervised cluster extraction through either 
  and constraint-based extraction (w/ a stability-based tiebreaker) or a mixed (weighted) 
  constraint and stability based objective extraction.
}
\usage{
extractFOSC(x, constraints = NA, alpha = 0, minPts = 2L, 
            prune_unstable = FALSE, 
            validate_constraints = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a valid hclust object.}
  \item{constraints}{ Either a list or matrix of pairwise constraints. If not supplied, an unsupervised measure of 'stability' is used to make local cuts and extract the 'optimal' clusters. See details. }
  \item{alpha}{ numeric; weight between [0, 1] for mixed-objective semisupervised extraction. Defaults to 0. See details. }
  \item{minPts}{ numeric; Defaults to 2. Only needed if class-less 'noise' is a valid 'label' in the model. }
  \item{prune_unstable}{ logical; should significantly unstable subtrees be pruned? The default
    is \code{FALSE} for the original optimal extraction framework (see Campello et al, 2013). See details for 
    what \code{TRUE} implies.}
  \item{validate_constraints}{ logical; should constraints be checked for validity? See details for what are 
  considered 'valid' constraints. }
}
\details{
Campello et al (2013) suggested a 'Framework for Optimal Selection of Clusters' (FOSC) as a framework to make local (non-horizontal) 
cuts to any cluster tree hierarchy. This function implements the original extraction algorithms as described by 
the framework for hclust objects. Traditional cluster extraction methods from hierarchical representations 
(such as 'hclust' objects) generally rely on global parameters or 'cutting' values which are used to partition a 
cluster hierarchy into a set of disjoint, 'flat' clusters. Although such methods are widespread, using global parameter 
settings are inherently limited in that they cannot capture patterns within the cluster hierarchy at varying \emph{local}
levels of granularity. 

Rather than partitioning a hierarchy based on the number of the cluster one expects to find (\eqn{k}) or based on some linkage 
distance threshold (\eqn{H}), the FOSC proposes that the most 'optimal' clusters may exist at varying distance thresholds in the
hierarchy. To enable this idea, FOSC requires one parameter (minPts) that represents \emph{the minimum number of 
points that constitute a valid cluster.} The first step of the FOSC algorithm is to traverse the given cluster hierarchy 
divisely, making recording 'new' clusters at each split if both branches are greater than or equal to minPts. Branches 
that contain less than minPts points at one or both branches inherit the parent clusters identity. Note that using FOSC, 
due to the restraint that minPts must be greater than or equal to 2, it is possible that the optimal cluster solution 
chosen makes local cuts that render parent branches of sizes less than minPts as 'noise' which are denoted as 0 in the 
final solution. 

Traversing the original cluster tree using minPts creates a new, 'simplified' cluster tree that is then post-processed 
recursively to extract clusters that maximize the following cost function for each cluster \eqn{C_i}{Ci}:

\deqn{\max_{\delta_2, \dots, \delta_k} J = \sum\limits_{i=2}^{k} \delta_i S(C_i)}{ J = \sum \delta S(Ci) for all i clusters }
where \eqn{S(C_i)}{S(Ci)} is the stability-based measure: 
\deqn{ S(C_i) = \sum_{x_j \in C_i}(\frac{1}{h_{min} (x_j, C_i)} - \frac{1}{h_{max} (C_i)}) }{ S(Ci) = \sum (1/Hmin(Xj, Ci) - 1/Hmax(Ci)) for all Xj in Ci}

and \eqn{\delta_i}{\delta} represents the indicator function, which constrains the solution space such that clusters must be 
disjoint (cannot assign more than 1 label to each cluster). The measure \eqn{S(C_i)}{S(Ci)} used by FOSC is an unsupervised 
validation measure based on the assumption that, if you vary the linkage/distance threshold across all possible values, 
more prominent clusters that survive over many threshold variations should be considered as stronger 
candidates of the 'optimal' solution. For this reason, using this measure to detect clusters is referred to as 
an unsupervised, \emph{stability-based} extraction approach. In some cases it may be useful to enact \emph{instance-level}
constraints that ensure the solution space conforms to linkage expectations known \emph{a priori}. This general idea 
of using preliminary expectations to augment the clustering solution will be referred to as \emph{semisupervised clustering}. 
If constraints are given in the call to 'extractFOSC', the following alternative objective function is maximized: 

\deqn{J = \frac{1}{2n_c}\sum\limits_{j=1}^n \gamma (x_j)}{J = 1/(2 * nc) \sum \gamma(Xj)}

where \eqn{n_c}{nc} is the total number of constraints given and \eqn{\gamma(x_j)}{\gamma(Xj)} represents the number of 
constraints involving object \eqn{x_j}{Xj} that are satisfied. In the case of ties (such as solutions where no constraints
were given), the unsupervised solution is used as a tiebreaker. See Campello et al (2013) for more details. 

As a third option, if one wishes to prioritize the degree at which the unsupervised and semisupervised solutions 
contribute to the overall optimal solution, the parameter \eqn{\alpha} can be set to enable the extraction of 
clusters that maxmize the following \emph{mixed} objective function: 

\deqn{J = \alpha S(C_i) + (1 - \alpha) \gamma(C_i))}{J = \alpha S(Ci) + (1 - \alpha) \gamma(Ci)}

FOSC expects the pairwise constraints to be passed in as either 1) an \eqn{n(n-1)/2} vector of integers representing 
the constraints, where 1 represents 'should-link', -1 represents 'should-not-link', and 0 represents 'no preference' 
(use unsupervised solution). See below for examples. Alternatively, if only a few constraints are needed, a named list 
representing the (symmetric!) adjacency list can be used, where the names correspond to indices of the points in the 
original data, and the values correspond to integer vectors of constraints (positive indices for 'should-link', negative 
indices for 'should-not-link'). Again, see the examples section for a demonstration of this.

The parameters to the input function correspond to the concepts discussed above. The \code{minPts} parameter to 
represent the minimum cluster size to extract. The optional \code{constraints} parameter contains the pairwise, 
instance-level constraints of the data. The optional \code{alpha} parameters controls whether the mixed objective 
function is used (if \code{alpha} is greater than 0). If the \code{validate_constraints} parameter is set to true, 
the constraints are checked (and fixed) for symmetry (if point A has a 'should-link' contraint with point B, point 
B should also have the same constraint). Asymmetric constraints are not supported. 

Unstable branch pruning was not discussed in the paper, however in some data sets it may be the case that specific 
subbranches scores are significantly greater than sibling and parent branches, and thus sibling branches should be 
considered as noise if their scores are cumulatively lower than the parents. This can happen in extremely nonhomogeneous 
data sets, where there exists locally very stable branches surrounded by 'unstable' but branches that however contain more than 
\emph{minPts} points. In such cases it may be useful to set this flag to true. 

}
\value{
  \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points (if any).}
  \item{hc }{The original hclust object augmented with the n-1 cluster-wide objective scores from the extraction 
  encoded in the 'stability', 'constraint', and 'total' named members.}
}
\references{
Campello, Ricardo JGB, Davoud Moulavi, Arthur Zimek, and Joerg Sander. "A framework for semi-supervised and unsupervised optimal extraction of clusters from hierarchies." \emph{Data Mining and Knowledge Discovery 27, no. 3 (2013): 344-371.}
}

\seealso{
\code{\link{hdbscan}}, \code{\link[stats]{cutree}}
}

\author{
  Matt Piekenbrock, 
  Michael Hahsler
}

\examples{
  data("moons")
  
  ## Regular HDBSCAN using stability-based extraction (unsupervised)
  cl <- hdbscan(moons, minPts = 5)
  print(cl$cluster)
  
  ## Constraint-based extraction from the HDBSCAN hierarchy 
  ## (w/ stability-based tiebreaker (semisupervised))
  cl_con <- extractFOSC(cl$hc, minPts = 5, constraints = list("12" = c(49, -47)))
  print(cl_con$cluster)
  
  ## Alternative formulation: Constraint-based extraction from the HDBSCAN hierarchy 
  ## (w/ stability-based tiebreaker (semisupervised)) using distance thresholds 
  dist_moons <- dist(moons)
  cl_con2 <- extractFOSC(cl$hc, minPts = 5, constraints = ifelse(dist_moons < 0.1, 1L, 
                                                          ifelse(dist_moons > 1, -1L, 0L)))
  print(cl_con2$cluster) # same as the second example 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ clustering }

