#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
#
# The __OPTICSXi__ R implementation was directly ported
# from the ELKI framework's (version 0.8.0) Java implementation with permission
# by the original author, Erich Schubert. This function is redistributed under
# the stricter [GNU AGPL v3](https://www.gnu.org/licenses/agpl-3.0.en.html).
#
# The package can be redistributed and modified under GNU AGPL v3.
# Remove the file and function to modify and redistribute the package
# under the GNU GPL v3 license.

#' @include optics.R
#' @rdname optics
#' @export
extractXi <-
  function(object,
           xi,
           minimum = FALSE,
           correctPredecessors = TRUE)
  {
    if (!inherits(object, "optics"))
      stop("extractXi only accepts xs resulting from dbscan::optics!")
    if (xi >= 1.0 ||
        xi <= 0.0)
      stop("The Xi parameter must be (0, 1)")

    SetOfClusters <-
      .xi_clusters(
        reach = object$reachdist[object$order],
        order = object$order,
        predecessor = object$predecessor,
        xi = xi,
        minPts = object$minPts,
        correctPredecessors = correctPredecessors
      )

    # Keep xi parameter, disable any previous flat clustering parameter
    object$xi <- xi
    object$eps_cl <- NA_real_

    # Zero-out clusters (only noise) if none found
    if (nrow(SetOfClusters) == 0) {
      warning(paste("No clusters were found with threshold:", xi))
      object$clusters_xi <- NULL
      object$cluster <- integer(length(object$order))
      return(invisible(object))
    }
    # Cluster data exists; organize it by starting and ending index, give arbitrary id
    object$clusters_xi <-
      data.frame(
        start = SetOfClusters[, 1L],
        end = SetOfClusters[, 2L],
        check.names = FALSE
      )
    object$clusters_xi <-
      object$clusters_xi[order(object$clusters_xi$start, object$clusters_xi$end), ]
    object$clusters_xi <-
      cbind(object$clusters_xi, list(cluster_id = seq_len(nrow(object$clusters_xi))))
    row.names(object$clusters_xi) <- NULL

    ## Populate cluster vector with either:
    ## 1. 'top-level' cluster labels to aid in plotting
    ## 2. 'local' or non-overlapping cluster labels if minimum == TRUE
    object$cluster <-
      extractClusterLabels(object$clusters_xi, object$order, minimum = minimum)

    # Remove non-local clusters if minimum was specified
    if (minimum) {
      object$clusters_xi <-
        object$clusters_xi[sort(unique(object$cluster))[-1], ]
    }

    class(object$cluster) <-
      unique(append(class(object$cluster), "xics"))
    class(object$clusters_xi) <-
      unique(append(class(object$clusters_xi), "xics"))
    object
  }

# Xi cluster extraction (Ankerst et al, 1999) with the filtering technique of
# Schubert and Gertz (2018).
#
# This is a port of OPTICSXi.extractClusters() in ELKI 0.8.0 (Java; 0-based
# indices there, 1-based here).
#
# reach: reachability distances in the cluster order
# order: the cluster order (point ids)
# predecessor: the predecessor of each point (indexed by point id)
#
# Returns a matrix with the columns start and end (positions in the order).
.xi_clusters <- function(reach, order, predecessor, xi, minPts, correctPredecessors = TRUE) {
  n <- length(order)
  ixi <- 1 - xi

  # steep up/down tests at position i for a steepness factor of ixi
  next_reach <- function(i) if (i < n) reach[i + 1L] else Inf
  steepUp <- function(i, ixi) {
    r <- reach[i]
    is.finite(r) && (i >= n || r <= reach[i + 1L] * ixi)
  }
  steepDown <- function(i, ixi) {
    nr <- next_reach(i)
    is.finite(nr) && nr <= reach[i] * ixi
  }

  # Remove obsolete steep down areas and update their maximum-in-between values
  update_sdaset <- function(mib, sdaset) {
    sdaset <- Filter(function(sda) sda$maximum * ixi > mib, sdaset)
    lapply(sdaset, function(sda) {
      if (mib > sda$mib)
        sda$mib <- mib
      sda
    })
  }

  # Remove trailing points from a cluster if their predecessor is not in the
  # cluster (Schubert and Gertz, 2018).
  predecessor_filter <- function(cstart, cend) {
    startval <- reach[cstart]
    while (cend > cstart) {
      # If the reachability is less than the one of the first point, then there
      # must be a point in between that "reached" the current point.
      if (reach[cend] < startval)
        break
      p <- predecessor[order[cend]]
      if (!is.na(p) && any(order[cstart:(cend - 1L)] == p))
        break
      cend <- cend - 1L
    }
    cend
  }

  mib <- 0
  sdaset <- list()
  clusters <- list()
  index <- 1L

  while (index <= n) {
    # Update the maximum-in-between
    mib <- max(mib, reach[index])

    # Steep down area
    if (steepDown(index, ixi)) {
      sdaset <- update_sdaset(mib, sdaset)
      startval <- reach[index]
      mib <- 0
      startsteep <- index
      endsteep <- index
      index <- index + 1L
      while (index <= n) {
        if (steepDown(index, ixi)) {
          endsteep <- index
        } else if (!steepDown(index, 1.0) || index - endsteep > minPts) {
          # not going downward at all or more than minPts flat steps
          break
        }
        index <- index + 1L
      }
      sdaset <- append(sdaset,
                       list(list(s = startsteep, e = endsteep,
                                 maximum = startval, mib = 0)))
      next
    }

    # Steep up area
    if (steepUp(index, ixi)) {
      sdaset <- update_sdaset(mib, sdaset)
      endsteep <- index
      mib <- reach[index]
      esuccr <- next_reach(index)
      # find the end of the steep up area, update mib again
      while (!is.infinite(esuccr) && index <= n) {
        index <- index + 1L
        if (index > n)
          break
        if (steepUp(index, ixi)) {
          endsteep <- index
          mib <- reach[index]
          esuccr <- next_reach(index)
        } else if (!steepUp(index, 1.0) || index - endsteep > minPts) {
          # not going upward at all or more than minPts flat steps
          break
        }
      }
      if (is.infinite(esuccr))
        index <- index + 1L

      # validate and compute clusters. Iterate backwards for the correct
      # hierarchy.
      for (sda in rev(sdaset)) {
        cstart <- sda$s
        cend <- endsteep
        if (correctPredecessors)
          cend <- predecessor_filter(cstart, cend)

        # end of the steep up area (use the next value)
        eU <- if (cend < n) reach[cend + 1L] else Inf

        # Condition 3b: maximum-in-between < end-of-steep-up * ixi
        if (sda$mib > min(sda$maximum, eU) * ixi)
          next

        # Condition 4 (case a is the default)
        if (sda$maximum * ixi >= eU) {
          # Case b
          while (cstart < cend && reach[cstart + 1L] * ixi > eU)
            cstart <- cstart + 1L
        } else if (eU * ixi >= sda$maximum) {
          # Case c
          while (cend > cstart && reach[cend] * ixi > sda$maximum)
            cend <- cend - 1L
        }

        # Filter again, mostly in case cstart was modified
        if (correctPredecessors)
          cend <- predecessor_filter(cstart, cend)

        # Condition 3a: obey minPts
        if (cend - cstart + 1L < minPts)
          next

        clusters <- append(clusters, list(c(as.numeric(cstart), as.numeric(cend))))
      }
      next
    }

    # flat: advance
    index <- index + 1L
  }

  if (length(clusters) == 0L)
    return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("start", "end"))))

  # remove duplicates
  res <- unique(do.call(rbind, clusters))
  colnames(res) <- c("start", "end")
  res
}
