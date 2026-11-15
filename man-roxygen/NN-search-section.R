# Shared description of the search parameters for kNN(), frNN() and sNN().
# Edit it here only. Use in a roxygen block with:
#   @templateVar type kNN    (or frNN / sNN)
#   @template NN-search-section
#
#' @details
#' ## Details on the search parameters:
#'
#' * `search` controls what method is used if `x` is a data matrix.
#'   Available methods are:
#'   - `"kdtree"` uses a Euclidean k-d tree to speed up NN search.
#'     This is the default with good overall performance.
#'   - `"bdtree"` uses a box-decomposition tree. Its spatial boxes can, in
#'     addition to splitting, also shrink. According to the ANN authors, this
#'     tree type may be faster than a k-d tree on tightly clustered data.
#'   - `"linear"` performs linear search directly on the data.
#'   - `"dist"` precomputes Euclidean distances first using [dist()]. While fast
#'     for small data sets, the memory requirements for the distance matrix grow
#'     with \eqn{O(n^2)}.
#'
#'   The trees and linear search use the highly optimized ANN library (Mount
#'   and Arya, 2010) for Euclidean distances.
#'
#'   If `x` is a distance matrix, then `search` is ignored and `"dist"`
#'   is used.
#'
#' * Tree specific parameters
#'   - `bucketSize`, `splitRule` and `shrinkRule` control how the tree is
#'     built. `bucketSize` and `splitRule` apply to k-d and bd trees;
#'     `shrinkRule` applies only to bd trees (`"none"` gives a tree that
#'     behaves like a k-d tree).
<% if (type == "frNN") { %>
#'   - `approx` uses the approximate fixed-radius search implemented in ANN.
#'     All neighbors up to a distance of `eps / (1 + approx)` are considered
#'     and none with a distance greater than `eps`. The points in between might
#'     be considered, so some actual neighbors may be omitted. This can lead
#'     to a significant speedup. The default `approx = 0` gives exact results.
#'     `approx` is only used by the tree-based searches. For more details see
#'     Mount and Arya (2010).
<% } else { %>
#'   - `approx` uses the approximate nearest neighbor search implemented in
#'     ANN. If `approx > 0`, the distance to the j-th returned neighbor is at
#'     most `(1 + approx)` times the distance to the true j-th nearest
#'     neighbor, so some true nearest neighbors may be omitted. This can lead
#'     to a significant speedup. The default `approx = 0` gives exact
#'     results. `approx` is only used by the tree-based searches. For more
#'     details see Mount and Arya (2010).
<% } %>
#'
#' ## Handling of `NA`s:
#'
#' The tree-based searches and `search = "linear"` cannot handle `NA`, `NaN` or
#' infinite values in the data (in `x` or `query`); an error is raised. (The ANN
#' library does not check its input and would silently return wrong results.)
#' With `search = "dist"`, the distances are calculated by [dist()], which
#' leaves out coordinates that are `NA` for a pair of points and scales up
#' the remaining sum. These distances are not exactly Euclidean distances.
#' The resulting distance matrix (and a [dist] object supplied as `x`)
#' cannot contain `NA`s, which happens, for example, if a point has only
#' `NA` values. Otherwise, remove or impute missing values before
#' performing the search.
