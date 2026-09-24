#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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

.ANNsplitRule <- c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST")

.matrixlike <- function(x) {
  if  (is.null(dim(x)))
       return(FALSE)

  # check that there is at least one row and one column!
  if (nrow(x) < 1L) stop("the provided data has 0 rows!")
  if (ncol(x) < 1L) stop("the provided data has 0 columns!")

  TRUE
}

.validate_integer_scalar <- function(x, name, min = 1) {
  valid <- length(x) == 1L &&
    is.numeric(x) &&
    typeof(x) %in% c("integer", "double") &&
    is.finite(x) &&
    x == trunc(x) &&
    x >= min

  if (!valid)
    stop(
      name,
      " must be a single, finite, integer-valued number >= ",
      min,
      ".",
      call. = FALSE
    )

  x
}

.validate_nonnegative_scalar <- function(x, name, allow_infinite = FALSE) {
  valid <- length(x) == 1L &&
    !is.na(x) &&
    is.numeric(x) &&
    typeof(x) %in% c("integer", "double") &&
    (is.finite(x) || (allow_infinite && x == Inf)) &&
    x >= 0

  if (!valid)
    stop(
      name,
      " must be a single, ",
      if (allow_infinite) "nonnegative" else "finite, nonnegative",
      " number.",
      call. = FALSE
    )

  x
}

.validate_bucket_size <- function(x) {
  .validate_integer_scalar(x, "bucketSize", min = 1)
}

.as_finite_numeric_matrix <- function(x, name = "x") {
  if (!.matrixlike(x))
    stop(name, " must be a matrix or data.frame.", call. = FALSE)

  x <- as.matrix(x)
  if (!is.numeric(x) || !typeof(x) %in% c("integer", "double"))
    stop(name, " must be a numeric matrix.", call. = FALSE)
  if (any(!is.finite(x)))
    stop(
      name,
      " cannot contain NA, NaN, or infinite values.",
      call. = FALSE
    )

  storage.mode(x) <- "double"
  x
}
