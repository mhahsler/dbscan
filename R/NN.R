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


adjacencylist <- function (x, ...) { UseMethod("adjacencylist", x) }
adjacencylist.frNN <- function(x, ...) x$id
adjacencylist.kNN <- function(x, ...)
  lapply(seq(nrow(x$id)), FUN = function(i) {
    ## filter NAs
    tmp <- x$id[i,]
    tmp[!is.na(tmp)]
    })

plot.NN <- function(x, data, main = NULL, ...) {
  if(is.null(main)) {
    if(inherits(x, "frNN")) main <- paste0("frNN graph (eps = ", x$eps, ")")
    if(inherits(x, "kNN")) main <- paste0(x$k, "-NN graph")
    if(inherits(x, "sNN")) main <- paste0("Shared NN graph (k=", x$k,
      ifelse(is.null(x$kt), "", paste0(", kt=", x$kt)), ")")
  }
  plot(data[,1:2], main = main, ...)

  id <- adjacencylist(x)

  ## use lines if it is from the same data
  ## FIXME: this test is not perfect, maybe we should have a parameter here or add the query points...
  if(length(id) == nrow(data)){
    for(i in 1:length(id)) {
      for(j in 1:length(id[[i]]))
        lines(x = c(data[i,1], data[id[[i]][j],1]),
          y = c(data[i,2], data[id[[i]][j],2]), ...)
    }
  } else { ## use colors if it was from a query
      for(i in 1:length(id)) {
    points(data[id[[i]],], col = i + 1L)
    }
  }
}

