library("dbscan")
library("testthat")

context("OPTICS")

# Use manual test to test against ELKI outputted files directly
manual_test <- FALSE
elki_out_location <- "~/Downloads/elki_xi_out/"
  
load(system.file("test_data/test_data.rda", package = "dbscan"))
load(system.file("test_data/elki_optics.rda", package = "dbscan"))

### run OPTICS XI xi=0.01, eps=1.0, minPts=5
x <- test_data
res <- optics(x, eps = 1.0,  minPts = 5)
res <- opticsXi(res, xi=0.01, minimum=F)

if (manual_test) { 
  ## Test against ELKI results manually 
  elki_optics_xi <- na.omit(data.frame(stringr::str_match(dir(elki_out_location), "cluster_(.*?)_(.*?).txt")[,2:3], stringsAsFactors=F))
  colnames(elki_optics_xi) <- c("start", "end")
  elki_optics_xi$start <- as.integer(elki_optics_xi$start)+1
  elki_optics_xi$end <- as.integer(elki_optics_xi$end)+1
  testthat::expect_identical(data.table::data.table(elki_optics_xi)[order(start, end)], res$clusters_xi[,list(start, end)][order(start, end)])
} else {
  load(system.file("test_data/elki_optics_xi.rda", package = "dbscan"))
  testthat::expect_identical(elki_optics_xi, data.frame(res$clusters_xi[,list(start, end)][order(start, end)]))
}

