#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2017 Michael Hahsler

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

#' Jarvis-Patrick Clustering
#'
#' Fast C++ implementation of the Jarvis-Patrick clustering which first builds
#' a shared nearest neighbor graph (k nearest neighbor sparsification) and then
#' places two points in the same cluster if they are in each other's nearest
#' neighbor list and they share at least kt nearest neighbors.
#'
#' Note: Following the original paper, the shared nearest neighbor list is
#' constructed as the k neighbors plus the point itself (as neighbor zero).
#' Therefore, the threshold \code{kt} can be in the range \eqn{[1, k]}.
#'
#' Fast nearest neighbors search with [kNN()] is only used if \code{x} is
#' a matrix. In this case Euclidean distance is used.
#'
#' @aliases jpclust print.general_clustering
#'
#' @param x a data matrix/data.frame (Euclidean distance is used), a
#' precomputed dist object or a kNN object created with [kNN()].
#' @param k Neighborhood size for nearest neighbor sparsification. If \code{x}
#' is a kNN object then \code{k} may be missing.
#' @param kt threshold on the number of shared nearest neighbors (including the
#' points themselves) to form clusters.
#' @param ...  additional arguments are passed on to the k nearest neighbor
#' search algorithm. See [kNN()] for details on how to control the
#' search strategy.
#'
#' @return A object of class `general_clustering` with the following
#' components:
#' \item{cluster }{A integer vector with cluster assignments. Zero
#' indicates noise points.}
#' \item{type }{ name of used clustering algorithm.}
#' \item{param }{ list of used clustering parameters. }
#'
#' @author Michael Hahsler
#' @seealso [dbscan()], [sNNclust()]
#' @references R. A. Jarvis and E. A. Patrick. 1973. Clustering Using a
#' Similarity Measure Based on Shared Near Neighbors. _IEEE Trans. Comput.
#' 22,_ 11 (November 1973), 1025-1034.
#' \doi{10.1109/T-C.1973.223640}
#' @keywords model clustering
#' @examples
#' data("DS3")
#'
#' # use a shared neighborhood of 20 points and require 12 shared neighbors
#' cl <- jpclust(DS3, k = 20, kt = 12)
#' cl
#'
#' plot(DS3, col = cl$cluster+1L, cex = .5)
#' # Note: JP clustering does not consider noise and thus,
#' # the sine wave points chain clusters together.
#'
#' # use a precomputed kNN object instead of the original data.
#' nn <- kNN(DS3, k = 30)
#' nn
#'
#' cl <- jpclust(nn, k = 20, kt = 12)
#' cl
#'
#' # cluster with noise removed (use low pointdensity to identify noise)
#' d <- pointdensity(DS3, eps = 25)
#' hist(d, breaks = 20)
#' DS3_noiseless <- DS3[d > 110,]
#'
#' cl <- jpclust(DS3_noiseless, k = 20, kt = 10)
#' cl
#'
#' plot(DS3_noiseless, col = cl$cluster+1L, cex = .5)
#'
#' @export jpclust
jpclust <- function(x, k, kt, ...) {
  # Create NN graph
  if (inherits(x, "kNN")) {
    if (missing(k))
      k <- nn$k
    nn <- x$id[, 1:k]
  } else {
    nn <- kNN(x, k, sort = FALSE, ...)$id
  }

  if (length(kt) != 1 || kt < 1 || kt > k)
    stop("kt needs to be a threshold in range [1, k].")

  # Perform clustering
  cl <- JP_int(nn, kt = as.integer(kt))

  structure(
    list(
      cluster = as.integer(factor(cl)),
      type = "Jarvis-Patrick clustering",
      param = list(k = k, kt = kt)
    ),
    class = c("general_clustering")
  )
}

print.general_clustering <- function(x, ...) {
  cl <- unique(x$cluster)
  cl <- length(cl[cl != 0L])

  writeLines(c(
    paste0(x$type, " for ", length(x$cluster), " objects."),
    paste0("Parameters: ",
      paste(
        names(x$param),
        unlist(x$param),
        sep = "=",
        collapse = ", "
      )),
    paste0(
      "The clustering contains ",
      cl,
      " cluster(s) and ",
      sum(x$cluster == 0L),
      " noise points."
    )
  ))

  print(table(x$cluster))
  cat("\n")

  writeLines(strwrap(paste0(
    "Available fields: ",
    paste(names(x), collapse = ", ")
  ), exdent = 18))
}
