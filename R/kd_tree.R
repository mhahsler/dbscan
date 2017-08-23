#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler, Matthew Piekenbrock

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

kd_tree <- function(x, bucketSize = 10L, splitRule = "suggest"){
  
  ## Check bucket size
  bucketSize <- if(missing(bucketSize)) 10L else as.integer(bucketSize)
  
  ## Check split rule 
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")
  
  ## Check data x 
  if(!.matrixlike(x)) stop("x needs to be a matrix")
  ## make sure x is numeric
  x <- as.matrix(x)
  if(storage.mode(x) == "integer") storage.mode(x) <- "double"
  if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")
  
  ## Get the ANN kd tree Xptr 
  kdtree_ptr <- kd_tree_int(x, bucketSize, splitRule)
  
  ## Sanity check
  if (deparse(kdtree_ptr) == "<pointer: 0x0>"){
    stop("Unable to create kd tree with the given data.")
  }
  
  ## Create kd tree R structure 
  res <- structure(list(bucketSize = bucketSize, 
                        splitRule = .ANNsplitRule[splitRule+1], 
                        call = match.call()), class = "kd_tree")
  attr(res, ".kd_tree_ptr") <- kdtree_ptr
  return(res)
}

print.kd_tree <- function(x){
  writeLines(c(
    paste0("KD Tree created from call: ", deparse(x$call)),
    paste0("Parameters: bucketSize = ", x$bucketSize, ", splitRule = ", x$splitRule)
  ))
}

summary.kd_tree <- function(x){
  kdtree_ptr <- attr(x, ".kd_tree_ptr")
  if (!is.null(kdtree_ptr) && class(kdtree_ptr) == "externalptr" && deparse(kdtree_ptr) != "<pointer: 0x0>"){
    printKdTreeStats(kdtree_ptr)
  }
}

str.kd_tree <- function(x){
  kdtree_ptr <- attr(x, ".kd_tree_ptr")
  if (!is.null(kdtree_ptr) && class(kdtree_ptr) == "externalptr" && deparse(kdtree_ptr) != "<pointer: 0x0>"){
    printKdTree(kdtree_ptr)
  }
}
