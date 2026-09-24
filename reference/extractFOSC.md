# Framework for the Optimal Extraction of Clusters from Hierarchies

Generic reimplementation of the *Framework for Optimal Selection of
Clusters* (FOSC; Campello et al, 2013) to extract clusterings from
hierarchical clustering (i.e.,
[hclust](https://rdrr.io/r/stats/hclust.html) objects). Can be
parameterized to perform unsupervised cluster extraction through a
stability-based measure, or semisupervised cluster extraction through
either a constraint-based extraction (with a stability-based tiebreaker)
or a mixed (weighted) constraint and stability-based objective
extraction.

## Usage

``` r
extractFOSC(
  x,
  constraints,
  alpha = 0,
  minPts = 2L,
  prune_unstable = FALSE,
  validate_constraints = FALSE
)
```

## Arguments

- x:

  a valid [hclust](https://rdrr.io/r/stats/hclust.html) object created
  via [`hclust()`](https://rdrr.io/r/stats/hclust.html) or
  [`hdbscan()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md).

- constraints:

  Either a list or matrix of pairwise constraints. If missing, an
  unsupervised measure of stability is used to make local cuts and
  extract the optimal clusters. See details.

- alpha:

  numeric; weight between \\\[0, 1\]\\ for mixed-objective
  semi-supervised extraction. Defaults to 0.

- minPts:

  numeric; Defaults to 2. Only needed if class-less noise is a valid
  label in the model.

- prune_unstable:

  logical; should significantly unstable subtrees be pruned? The default
  is `FALSE` for the original optimal extraction framework (see Campello
  et al, 2013). See details for what `TRUE` implies.

- validate_constraints:

  logical; should constraints be checked for validity? See details for
  what are considered valid constraints.

## Value

A list with the elements:

- cluster :

  A integer vector with cluster assignments. Zero indicates noise points
  (if any).

- hc :

  The original [hclust](https://rdrr.io/r/stats/hclust.html) object with
  additional list elements `"stability"`, `"constraint"`, and `"total"`
  for the \\n - 1\\ cluster-wide objective scores from the extraction.

## Details

Campello et al (2013) suggested a *Framework for Optimal Selection of
Clusters* (FOSC) as a framework to make local (non-horizontal) cuts to
any cluster tree hierarchy. This function implements the original
extraction algorithms as described by the framework for hclust objects.
Traditional cluster extraction methods from hierarchical representations
(such as [hclust](https://rdrr.io/r/stats/hclust.html) objects)
generally rely on global parameters or cutting values which are used to
partition a cluster hierarchy into a set of disjoint, flat clusters.
This is implemented in R in function
[`stats::cutree()`](https://rdrr.io/r/stats/cutree.html). Although such
methods are widespread, using global parameter settings are inherently
limited in that they cannot capture patterns within the cluster
hierarchy at varying *local* levels of granularity.

Rather than partitioning a hierarchy based on the number of the cluster
one expects to find (\\k\\) or based on some linkage distance threshold
(\\H\\), the FOSC proposes that the optimal clusters may exist at
varying distance thresholds in the hierarchy. To enable this idea, FOSC
requires one parameter (minPts) that represents *the minimum number of
points that constitute a valid cluster.* The first step of the FOSC
algorithm is to traverse the given cluster hierarchy divisively,
recording new clusters at each split if both branches represent more
than or equal to minPts. Branches that contain less than minPts points
at one or both branches inherit the parent clusters identity. Note that
using FOSC, due to the constraint that minPts must be greater than or
equal to 2, it is possible that the optimal cluster solution chosen
makes local cuts that render parent branches of sizes less than minPts
as noise, which are denoted as 0 in the final solution.

Traversing the original cluster tree using minPts creates a new,
simplified cluster tree that is then post-processed recursively to
extract clusters that maximize for each cluster \\C_i\\ the cost
function

\$\$\max\_{\delta_2, \dots, \delta_k} J = \sum\limits\_{i=2}^{k}
\delta_i S(C_i)\$\$ where \\S(C_i)\\ is the stability-based measure as
\$\$ S(C_i) = \sum\_{x_j \in C_i}(\frac{1}{h\_{min} (x_j, C_i)} -
\frac{1}{h\_{max} (C_i)}) \$\$

\\\delta_i\\ represents an indicator function, which constrains the
solution space such that clusters must be disjoint (cannot assign more
than 1 label to each cluster). The measure \\S(C_i)\\ used by FOSC is an
unsupervised validation measure based on the assumption that, if you
vary the linkage/distance threshold across all possible values, more
prominent clusters that survive over many threshold variations should be
considered as stronger candidates of the optimal solution. For this
reason, using this measure to detect clusters is referred to as an
unsupervised, *stability-based* extraction approach. In some cases it
may be useful to enact *instance-level* constraints that ensure the
solution space conforms to linkage expectations known *a priori*. This
general idea of using preliminary expectations to augment the clustering
solution will be referred to as *semisupervised clustering*. If
constraints are given in the call to `extractFOSC()`, the following
alternative objective function is maximized:

\$\$J = \frac{1}{2n_c}\sum\limits\_{j=1}^n \gamma (x_j)\$\$

\\n_c\\ is the total number of constraints given and \\\gamma(x_j)\\
represents the number of constraints involving object \\x_j\\ that are
satisfied. In the case of ties (such as solutions where no constraints
were given), the unsupervised solution is used as a tiebreaker. See
Campello et al (2013) for more details.

As a third option, if one wishes to prioritize the degree at which the
unsupervised and semisupervised solutions contribute to the overall
optimal solution, the parameter \\\alpha\\ can be set to enable the
extraction of clusters that maximize the `mixed` objective function

\$\$J = \alpha S(C_i) + (1 - \alpha) \gamma(C_i))\$\$

FOSC expects the pairwise constraints to be passed as either 1) an
\\n(n-1)/2\\ vector of integers representing the constraints, where 1
represents should-link, -1 represents should-not-link, and 0 represents
no preference using the unsupervised solution (see below for examples).
Alternatively, if only a few constraints are needed, a named list
representing the (symmetric) adjacency list can be used, where the names
correspond to indices of the points in the original data, and the values
correspond to integer vectors of constraints (positive indices for
should-link, negative indices for should-not-link). Again, see the
examples section for a demonstration of this.

The parameters to the input function correspond to the concepts
discussed above. The `minPts` parameter to represent the minimum cluster
size to extract. The optional `constraints` parameter contains the
pairwise, instance-level constraints of the data. The optional `alpha`
parameters controls whether the mixed objective function is used (if
`alpha` is greater than 0). If the `validate_constraints` parameter is
set to true, the constraints are checked (and fixed) for symmetry (if
point A has a should-link constraint with point B, point B should also
have the same constraint). Asymmetric constraints are not supported.

Unstable branch pruning was not discussed by Campello et al (2013),
however in some data sets it may be the case that specific subbranches
scores are significantly greater than sibling and parent branches, and
thus sibling branches should be considered as noise if their scores are
cumulatively lower than the parents. This can happen in extremely
nonhomogeneous data sets, where there exists locally very stable
branches surrounded by unstable branches that contain more than `minPts`
points. `prune_unstable = TRUE` will remove the unstable branches.

## References

Campello, Ricardo JGB, Davoud Moulavi, Arthur Zimek, and Joerg Sander
(2013). A framework for semi-supervised and unsupervised optimal
extraction of clusters from hierarchies. *Data Mining and Knowledge
Discovery* 27(3): 344-371.
[doi:10.1007/s10618-013-0311-4](https://doi.org/10.1007/s10618-013-0311-4)

## See also

[`hclust()`](https://rdrr.io/r/stats/hclust.html),
[`hdbscan()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md),
[`stats::cutree()`](https://rdrr.io/r/stats/cutree.html)

Other clustering functions:
[`dbscan()`](http://michael.hahsler.net/dbscan/reference/dbscan.md),
[`hdbscan()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md),
[`jpclust()`](http://michael.hahsler.net/dbscan/reference/jpclust.md),
[`ncluster()`](http://michael.hahsler.net/dbscan/reference/ncluster.md),
[`optics()`](http://michael.hahsler.net/dbscan/reference/optics.md),
[`sNNclust()`](http://michael.hahsler.net/dbscan/reference/sNNclust.md)

## Author

Matt Piekenbrock

## Examples

``` r
data("moons")

## Regular HDBSCAN using stability-based extraction (unsupervised)
cl <- hdbscan(moons, minPts = 5)
cl$cluster
#>   [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
#>  [38] 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 1 1 2 2 1 2 2 2 1 1 2 1 1 1 2 1 1 2 2
#>  [75] 2 1 2 1 2 2 1 2 2 1 1 2 2 1 2 2 1 2 2 2 1 1 1 2 1 1

## Constraint-based extraction from the HDBSCAN hierarchy
## (w/ stability-based tiebreaker (semisupervised))
cl_con <- extractFOSC(cl$hc, minPts = 5,
  constraints = list("12" = c(49, -47)))
cl_con$cluster
#>   [1] 6 5 5 5 5 6 5 6 6 5 6 0 6 6 5 5 5 5 5 5 5 6 5 6 5 5 6 6 6 5 6 0 5 6 6 6 5
#>  [38] 5 0 6 6 6 5 6 5 6 5 6 0 5 3 3 1 1 1 1 1 3 3 1 3 3 3 1 1 3 1 1 1 3 1 1 3 3
#>  [75] 3 1 3 1 3 3 1 3 3 1 1 3 3 1 3 3 1 3 3 3 1 1 1 3 1 1

## Alternative formulation: Constraint-based extraction from the HDBSCAN hierarchy
## (w/ stability-based tiebreaker (semisupervised)) using distance thresholds
dist_moons <- dist(moons)
cl_con2 <- extractFOSC(cl$hc, minPts = 5,
  constraints = ifelse(dist_moons < 0.1, 1L,
                ifelse(dist_moons > 1, -1L, 0L)))

cl_con2$cluster # same as the second example
#>   [1] 6 5 5 5 5 6 5 6 6 5 6 0 6 6 5 5 5 5 5 5 5 6 5 6 5 5 6 6 6 5 6 0 5 6 6 6 5
#>  [38] 5 0 6 6 6 5 6 5 6 5 6 0 5 3 3 1 1 1 1 1 3 3 1 3 3 3 1 1 3 1 1 1 3 1 1 3 3
#>  [75] 3 1 3 1 3 3 1 3 3 1 1 3 3 1 3 3 1 3 3 3 1 1 1 3 1 1
```
