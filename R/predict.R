.predict_frNN <- function(newdata, data, clusters, eps, ...) {
  if (is.null(newdata)) return(clusters)

  nn <- frNN(rbind(data, newdata), eps = eps,
    sort = TRUE, ...)$id[-(1:nrow(data))]

  sapply(nn, function(x) {
    x <- x[x<=nrow(data)]
    x <- clusters[x][x>0][1]
    x[is.na(x)] <- 0L
    x
  })
}

predict.dbscan_fast <- function (object, newdata = NULL, data, ...)
  .predict_frNN(newdata, data, object$cluster, object$eps, ...)

predict.optics <- function (object, newdata = NULL, data, ...) {
  if (is.null(object$cluster)) stop("no extracted clustering available in object! run extractDBSCAN or extractXi first.")
  .predict_frNN(newdata, data, object$cluster, object$eps_cl, ...)
}
