# Shared parameter documentation for the neighbor search functions kNN(),
# frNN() and sNN(). Edit it here only. Use in a roxygen block with:
#   @templateVar type kNN    (or frNN / sNN)
#   @template NN-search-params
#
#' @param search nearest neighbor search strategy (one of `"kdtree"`, `"linear"`,
#' `"dist"` or `"bdtree"`). See the section on the search parameters.
#' @param bucketSize max size of the k-d or bd tree leafs.
#' @param splitRule rule to split the k-d or bd tree. One of `"std"`, `"midpt"`, `"fair"`,
#' `"sl_midpt"`, `"sl_fair"` or `"suggest"` (SL stands for sliding). The default
#' value `"suggest"` uses ANNs best guess.
#' @param shrinkRule shrinking rule used to build the box-decomposition tree
#' (only used for `search = "bdtree"`). One of `"none"` (no shrinking, equivalent
#' to a k-d tree), `"simple"`, `"centroid"` or `"suggest"`. The default value
#' `"suggest"` uses ANNs best guess.
<% if (type == "frNN") { %>
#' @param approx a non-negative number used for approximate nearest neighbors.
#' All NN up to a distance of `eps / (1 + approx)` are considered and none with
#' a distance greater than `eps`. The points in between might be considered.
#' Some actual NN may be omitted, but the search will be faster.
#' Ignored by `search = "linear"` and `"dist"`.
<% } else { %>
#' @param approx approximation factor for the tree-based searches. The j-th
#' neighbor found is at most `(1 + approx)` times as far away as the true j-th
#' nearest neighbor. `0` (default) gives exact results; larger values are
#' faster, but some true nearest neighbors may be omitted. Ignored by
#' `search = "linear"` and `"dist"`.
<% } %>
