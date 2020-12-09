# ANN uses a global KD_TRIVIAL structure which needs to be removed.
.onUnload <- function(libpath) {
  ANN_cleanup()
  #cat("Cleaning up after ANN.\n")
}
