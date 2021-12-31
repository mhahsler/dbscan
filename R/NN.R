#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

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

#' NN --- Nearest Neighbors Superclass
#'
#' NN is an abstract S3 superclass for the classes of the objects returned
#' by [kNN()], [frNN()] and [sNN()]. Methods for sorting, plotting and getting an
#' adjacency list are defined.
#'
#' @name NN
#' @aliases NN
#' @family NN functions
#'
#' @param x a `NN` object
#' @param ... further parameters
#' @param decreasing sort in decreasing order?
#' @param data that was used to create `x`
#' @param main title
#'
#' @section Subclasses:
#' [kNN], [frNN] and [sNN]
#'
#' @author Michael Hahsler
#' @keywords model
#' @examples
#' data(iris)
#' x <- iris[, -5]
#'
#' # finding kNN directly in data (using a kd-tree)
#' nn <- kNN(x, k=5)
#' nn
#'
#' # plot the kNN where NN are shown as line conecting points.
#' plot(nn, x)
#'
#' # show the first few elements of the adjacency list
#' head(adjacencylist(nn))
#'
#' # create a graph and find connected components (if igraph is installed)
#' if("igraph" %in% installed.packages()){
#'   library("igraph")
#'   g <- graph_from_adj_list(adjacencylist(nn))
#'   comp <- components(g)
#'   plot(x, col = comp$membership)
#'
#'   # detect clusters (communities) with the label propagation algorithm
#'   cl <- membership(cluster_label_prop(g))
#'   plot(x, col = cl)
#' }
NULL

#' @rdname NN
adjacencylist <- function (x, ...)
  UseMethod("adjacencylist", x)

#' @rdname NN
adjacencylist.NN <- function (x, ...) {
  stop("needs to be implemented by a subclass")
  }

#' @rdname NN
sort.NN <- function(x, decreasing = FALSE, ...) {
  stop("needs to be implemented by a subclass")
  }


#' @rdname NN
plot.NN <- function(x, data, main = NULL, ...) {
  if (is.null(main)) {
    if (inherits(x, "frNN"))
      main <- paste0("frNN graph (eps = ", x$eps, ")")
    if (inherits(x, "kNN"))
      main <- paste0(x$k, "-NN graph")
    if (inherits(x, "sNN"))
      main <- paste0("Shared NN graph (k=", x$k,
        ifelse(is.null(x$kt), "", paste0(", kt=", x$kt)), ")")
  }
  plot(data[, 1:2], main = main, ...)

  id <- adjacencylist(x)

  ## use lines if it is from the same data
  ## FIXME: this test is not perfect, maybe we should have a parameter here or add the query points...
  if (length(id) == nrow(data)) {
    for (i in 1:length(id)) {
      for (j in 1:length(id[[i]]))
        lines(x = c(data[i, 1], data[id[[i]][j], 1]),
          y = c(data[i, 2], data[id[[i]][j], 2]),
          ...)
    }
  } else {
    ## use colors if it was from a query
    for (i in 1:length(id)) {
      points(data[id[[i]], ], col = i + 1L)
    }
  }
}
