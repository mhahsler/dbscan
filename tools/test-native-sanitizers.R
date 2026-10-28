args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(args[1], .libPaths()))
library(dbscan)
stopifnot(normalizePath(find.package("dbscan")) ==
            normalizePath(file.path(args[1], "dbscan")))

# Deliver real SIGINTs during each ANN routine, rather than relying on timing
# assertions in the ordinary tests. Linear searches keep these calls running
# long enough for the child to signal the parent. Each child is reaped before
# proceeding, and subsequent calls exercise the package after unwinding.
set.seed(20261029)
x <- matrix(runif(90000), ncol = 3)
parent <- Sys.getpid()
operations <- list(
  function() kNN(x, 5, search = "linear"),
  function() kNN(x, 5, query = x, search = "linear"),
  function() frNN(x, .001, search = "linear"),
  function() frNN(x, .001, query = x, search = "linear"),
  function() dbscan(x, .001, 3, search = "linear"),
  function() lof(x, 5, search = "linear"),
  function() optics(x, .001, 3, search = "linear"),
  function() pointdensity(x, .001, search = "linear")
)
for (operation in operations) {
  sender <- parallel::mcparallel({
    Sys.sleep(.1)
    tools::pskill(parent, signal = 2L)
    NULL
  }, mc.set.seed = FALSE)
  interrupted <- tryCatch({ operation(); FALSE }, interrupt = function(e) TRUE,
                          finally = parallel::mccollect(sender))
  stopifnot(interrupted)
}
rm(x, operations, operation, sender)
cat("Native interrupt checks passed.\n")

testthat::test_dir(
  args[2],
  package = "dbscan",
  load_package = "installed",
  reporter = "check",
  stop_on_failure = TRUE
)

# Release ANN's package-global KD_TRIVIAL before the exit-time leak scan.
unloadNamespace("dbscan")
invisible(gc())
