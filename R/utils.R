# Base R introduced the %||% null-coalescing operator in R version 4.4.0
# Backport for 4.3.0

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
