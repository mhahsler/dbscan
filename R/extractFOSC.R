#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler, Matt Piekenbrock

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#' Framework for Optimal Selection of Clusters
#'
#' Generic reimplementation of the Framework for Optimal Selection of Clusters
#' (FOSC; Campello et al, 2013). Can be parameterized to perform unsupervised
#' cluster extraction through a stability-based measure, or semisupervised
#' cluster extraction through either a constraint-based extraction (with a
#' stability-based tiebreaker) or a mixed (weighted) constraint and
#' stability-based objective extraction.
#'
#' Campello et al (2013) suggested a 'Framework for Optimal Selection of
#' Clusters' (FOSC) as a framework to make local (non-horizontal) cuts to any
#' cluster tree hierarchy. This function implements the original extraction
#' algorithms as described by the framework for hclust objects. Traditional
#' cluster extraction methods from hierarchical representations (such as
#' 'hclust' objects) generally rely on global parameters or cutting values
#' which are used to partition a cluster hierarchy into a set of disjoint, flat
#' clusters. Although such methods are widespread, using global parameter
#' settings are inherently limited in that they cannot capture patterns within
#' the cluster hierarchy at varying _local_ levels of granularity.
#'
#' Rather than partitioning a hierarchy based on the number of the cluster one
#' expects to find (\eqn{k}) or based on some linkage distance threshold
#' (\eqn{H}), the FOSC proposes that the optimal clusters may exist at varying
#' distance thresholds in the hierarchy. To enable this idea, FOSC requires one
#' parameter (minPts) that represents _the minimum number of points that
#' constitute a valid cluster._ The first step of the FOSC algorithm is to
#' traverse the given cluster hierarchy divisively, recording new clusters at
#' each split if both branches represent more than or equal to minPts. Branches
#' that contain less than minPts points at one or both branches inherit the
#' parent clusters identity. Note that using FOSC, due to the constraint that
#' minPts must be greater than or equal to 2, it is possible that the optimal
#' cluster solution chosen makes local cuts that render parent branches of
#' sizes less than minPts as noise, which are denoted as 0 in the final
#' solution.
#'
#' Traversing the original cluster tree using minPts creates a new, simplified
#' cluster tree that is then post-processed recursively to extract clusters
#' that maximize for each cluster \eqn{C_i}{Ci} the cost function
#'
#' \deqn{\max_{\delta_2, \dots, \delta_k} J = \sum\limits_{i=2}^{k} \delta_i
#' S(C_i)}{ J = \sum \delta S(Ci) for all i clusters, } where
#' \eqn{S(C_i)}{S(Ci)} is the stability-based measure as \deqn{ S(C_i) =
#' \sum_{x_j \in C_i}(\frac{1}{h_{min} (x_j, C_i)} - \frac{1}{h_{max} (C_i)})
#' }{ S(Ci) = \sum (1/Hmin(Xj, Ci) - 1/Hmax(Ci)) for all Xj in Ci.}
#'
#' \eqn{\delta_i}{\delta} represents an indicator function, which constrains
#' the solution space such that clusters must be disjoint (cannot assign more
#' than 1 label to each cluster). The measure \eqn{S(C_i)}{S(Ci)} used by FOSC
#' is an unsupervised validation measure based on the assumption that, if you
#' vary the linkage/distance threshold across all possible values, more
#' prominent clusters that survive over many threshold variations should be
#' considered as stronger candidates of the optimal solution. For this reason,
#' using this measure to detect clusters is referred to as an unsupervised,
#' _stability-based_ extraction approach. In some cases it may be useful
#' to enact _instance-level_ constraints that ensure the solution space
#' conforms to linkage expectations known _a priori_. This general idea of
#' using preliminary expectations to augment the clustering solution will be
#' referred to as _semisupervised clustering_. If constraints are given in
#' the call to `extractFOSC()`, the following alternative objective function
#' is maximized:
#'
#' \deqn{J = \frac{1}{2n_c}\sum\limits_{j=1}^n \gamma (x_j)}{J = 1/(2 * nc)
#' \sum \gamma(Xj)}
#'
#' \eqn{n_c}{nc} is the total number of constraints given and
#' \eqn{\gamma(x_j)}{\gamma(Xj)} represents the number of constraints involving
#' object \eqn{x_j}{Xj} that are satisfied. In the case of ties (such as
#' solutions where no constraints were given), the unsupervised solution is
#' used as a tiebreaker. See Campello et al (2013) for more details.
#'
#' As a third option, if one wishes to prioritize the degree at which the
#' unsupervised and semisupervised solutions contribute to the overall optimal
#' solution, the parameter \eqn{\alpha} can be set to enable the extraction of
#' clusters that maximize the `mixed` objective function
#'
#' \deqn{J = \alpha S(C_i) + (1 - \alpha) \gamma(C_i))}{J = \alpha S(Ci) + (1 -
#' \alpha) \gamma(Ci).}
#'
#' FOSC expects the pairwise constraints to be passed as either 1) an
#' \eqn{n(n-1)/2} vector of integers representing the constraints, where 1
#' represents should-link, -1 represents should-not-link, and 0 represents no
#' preference using the unsupervised solution (see below for examples).
#' Alternatively, if only a few constraints are needed, a named list
#' representing the (symmetric) adjacency list can be used, where the names
#' correspond to indices of the points in the original data, and the values
#' correspond to integer vectors of constraints (positive indices for
#' should-link, negative indices for should-not-link). Again, see the examples
#' section for a demonstration of this.
#'
#' The parameters to the input function correspond to the concepts discussed
#' above. The \code{minPts} parameter to represent the minimum cluster size to
#' extract. The optional \code{constraints} parameter contains the pairwise,
#' instance-level constraints of the data. The optional \code{alpha} parameters
#' controls whether the mixed objective function is used (if \code{alpha} is
#' greater than 0). If the \code{validate_constraints} parameter is set to
#' true, the constraints are checked (and fixed) for symmetry (if point A has a
#' should-link constraint with point B, point B should also have the same
#' constraint). Asymmetric constraints are not supported.
#'
#' Unstable branch pruning was not discussed by Campello et al (2013), however
#' in some data sets it may be the case that specific subbranches scores are
#' significantly greater than sibling and parent branches, and thus sibling
#' branches should be considered as noise if their scores are cumulatively
#' lower than the parents. This can happen in extremely nonhomogeneous data
#' sets, where there exists locally very stable branches surrounded by unstable
#' branches that contain more than `minPts` points. \code{prune_unstable =
#' TRUE} will remove the unstable branches.
#'
#' @param x a valid [hclust] object.
#' @param constraints Either a list or matrix of pairwise constraints. If
#' missing, an unsupervised measure of stability is used to make local cuts and
#' extract the optimal clusters. See details.
#' @param alpha numeric; weight between \eqn{[0, 1]} for mixed-objective
#' semi-supervised extraction. Defaults to 0.
#' @param minPts numeric; Defaults to 2. Only needed if class-less noise is a
#' valid label in the model.
#' @param prune_unstable logical; should significantly unstable subtrees be
#' pruned? The default is \code{FALSE} for the original optimal extraction
#' framework (see Campello et al, 2013). See details for what \code{TRUE}
#' implies.
#' @param validate_constraints logical; should constraints be checked for
#' validity? See details for what are considered valid constraints.
#' @return A list with the elements:
#'
#' \item{cluster }{A integer vector with cluster assignments. Zero
#' indicates noise points (if any).}
#' \item{hc }{The original [hclust] object
#' augmented with the n-1 cluster-wide objective scores from the extraction
#' encoded in the 'stability', 'constraint', and 'total' named members.}
#'
#' @author Matt Piekenbrock
#' @seealso [hdbscan()], [stats::cutree()]
#' @references Campello, Ricardo JGB, Davoud Moulavi, Arthur Zimek, and Joerg
#' Sander (2013). "A framework for semi-supervised and unsupervised optimal
#' extraction of clusters from hierarchies." _Data Mining and Knowledge
#' Discovery_ 27(3): 344-371.
#' \doi{10.1007/s10618-013-0311-4}
#' @keywords model clustering
#' @examples
#' data("moons")
#'
#' ## Regular HDBSCAN using stability-based extraction (unsupervised)
#' cl <- hdbscan(moons, minPts = 5)
#' cl$cluster
#'
#' ## Constraint-based extraction from the HDBSCAN hierarchy
#' ## (w/ stability-based tiebreaker (semisupervised))
#' cl_con <- extractFOSC(cl$hc, minPts = 5,
#'   constraints = list("12" = c(49, -47)))
#' cl_con$cluster
#'
#' ## Alternative formulation: Constraint-based extraction from the HDBSCAN hierarchy
#' ## (w/ stability-based tiebreaker (semisupervised)) using distance thresholds
#' dist_moons <- dist(moons)
#' cl_con2 <- extractFOSC(cl$hc, minPts = 5,
#'   constraints = ifelse(dist_moons < 0.1, 1L,
#'                 ifelse(dist_moons > 1, -1L, 0L)))
#'
#' cl_con2$cluster # same as the second example
#'
#' @export extractFOSC
extractFOSC <-
  function(x,
    constraints,
    alpha = 0,
    minPts = 2L,
    prune_unstable = FALSE,
    validate_constraints = FALSE) {
    if (!inherits(x, "hclust"))
      stop("extractFOSC expects 'x' to be a valid hclust object.")

    # if constraints are given then they need to be a list, a matrix or a vector
    if (!(
      missing(constraints) ||
        is.list(constraints) ||
        is.matrix(constraints) ||
        is.numeric(constraints)
    ))
      stop("extractFOSC expects constraints to be either an adjacency list or adjacency matrix.")

    if (!minPts >= 2)
      stop("minPts must be at least 2.")
    if (alpha < 0 ||
        alpha > 1)
      stop("alpha can only takes values between [0, 1].")
    n <- nrow(x$merge) + 1L

    ## First step for both unsupervised and semisupervised - compute stability scores
    cl_tree <- computeStability(x, minPts)

    ## Unsupervised Extraction
    if (missing(constraints)) {
      cl_tree <- extractUnsupervised(cl_tree, prune_unstable)
    }
    ## Semi-supervised Extraction
    else {
      ## If given as adjacency-list form
      if (is.list(constraints)) {
        ## Checks for proper indexing, symmetry of constraints, etc.
        if (validate_constraints) {
          is_valid <- max(as.integer(names(constraints))) < n
          is_valid <-
            is_valid && all(sapply(constraints, function(ilc)
              all(ilc <= n)))
          if (!is_valid) {
            stop("Detected constraint indices not in the interval [1, n]")
          }
          constraints <- validateConstraintList(constraints, n)
        }
        cl_tree <-
          extractSemiSupervised(cl_tree, constraints, alpha, prune_unstable)
      }
      ## Adjacency matrix given (probably from dist object), retrieve adjacency list form
      else if (is.vector(constraints)) {
        if (!all(constraints %in% c(-1, 0, 1))) {
          stop(
            "'extractFOSC' only accepts instance-level constraints. See ?extractFOSC for more details."
          )
        }
        ## Checks for proper integer labels, symmetry of constraints, length of vector, etc.
        if (validate_constraints) {
          is_valid <- length(constraints) == choose(n, 2)
          constraints_list <-
            validateConstraintList(distToAdjacency(constraints, n), n)
        } else {
          constraints_list <-  distToAdjacency(constraints, n)
        }
        cl_tree <-
          extractSemiSupervised(cl_tree, constraints_list, alpha, prune_unstable)
      }
      ## Full nxn adjacency-matrix given, give warning and retrieve adjacency list form
      else if (is.matrix(constraints)) {
        if (!all(constraints %in% c(-1, 0, 1))) {
          stop(
            "'extractFOSC' only accepts instance-level constraints. See ?extractFOSC for more details."
          )
        }
        if (!all(dim(constraints) == c(n, n))) {
          stop("Given matrix is not square.")
        }
        warning(
          "Full nxn matrix given; extractFOCS does not support asymmetric relational constraints. Using lower triangular."
        )

        constraints <- constraints[lower.tri(constraints)]

        ## Checks for proper integer labels, symmetry of constraints, length of vector, etc.
        if (validate_constraints) {
          is_valid <- length(constraints) == choose(n, 2)
          constraints_list <-
            validateConstraintList(distToAdjacency(constraints, n), n)
        } else {
          constraints_list <- distToAdjacency(constraints, n)
        }
        cl_tree <-
          extractSemiSupervised(cl_tree, constraints_list, alpha, prune_unstable)
      } else {
        stop(paste(
          "'extractFOSC' doesn't know how to handle constraints of type",
          class(constraints)
        ))
      }
    }
    total_stab <-
      if (is.null(attr(cl_tree, "total_stability")))
        1
    else
      attr(cl_tree, "total_stability")
    cl_track <- attr(cl_tree, "cl_tracker")
    stability_score <-
      unlist(sapply(cl_track, function(cid)
        cl_tree[[as.character(cid)]]$stability))
    constraint_score <-
      unlist(sapply(cl_track, function(cid)
        cl_tree[[as.character(cid)]]$vscore))
    total_score <-
      unlist(sapply(cl_track, function(cid)
        cl_tree[[as.character(cid)]]$vscore))
    out <- append(
      x,
      list(
        "cluster" = cl_track,
        "stability" = stability_score,
        "constraint" = constraint_score,
        "total" = total_score
      )
    )
    extraction_type <-
      ifelse(
        missing(constraints),
        "(w/ stability-based extraction)",
        ifelse(
          alpha == 0,
          "(w/ constraint-based extraction)",
          "(w/ mixed-objective extraction)"
        )
      )
    substrs <- unlist(strsplit(x$method, split = " \\(w\\/"))
    out[["method"]] <-
      if (length(substrs) > 1)
        paste(substrs[[1]], extraction_type)
    else
      paste(out[["method"]], extraction_type)
    class(out) <- "hclust"
    return(list(cluster = attr(cl_tree, "cluster"), hc = out))
  }
