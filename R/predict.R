.predict_frNN <- function(newdata, data, clusters, eps, ...) {
  if (is.null(newdata)) return(clusters)

  # calculate the frNN between newdata and data (only keep entries for newdata)
  nn <- frNN(rbind(data, newdata), eps = eps,
    sort = TRUE, ...)$id[-(1:nrow(data))]

  # get rid of indices of newdata and pick first neighbor
  cl <- sapply(nn, function(x) clusters[x[x<=nrow(data) & x>0][1]])
  cl[is.na(cl)] <- 0L
  cl
}

predict.dbscan_fast <- function (object, newdata = NULL, data, ...)
  .predict_frNN(newdata, data, object$cluster, object$eps, ...)

predict.optics <- function (object, newdata = NULL, data, ...) {
  if (is.null(object$cluster)) stop("no extracted clustering available in object! run extractDBSCAN or extractXi first.")
  .predict_frNN(newdata, data, object$cluster, object$eps_cl, ...)
}

predict.hdbscan <- function(object, newdata = NULL, data, ...) {
  
  k = 2 * object$minPts
  
  # get all the nearest neighbor IDs for newdata and their distances
  nn <- kNN(rbind(data, newdata), k = k, sort = TRUE)$id[-(1:nrow(data)),]
  nn_dist <- kNNdist(rbind(data, newdata), k = k, all = TRUE)[-(1:nrow(data)),]
  core_dist <- nn_dist[,object$minPts - 1]
  
  # for each new data point, compute MRD for its nearest neighbors
  nn_mrd <- mrd_m(nn_dist, core_dist)
  
  # reorder the neighbors by their MRD
  nn_order <- t(apply(nn_mrd, 1, order))
  
  # get the cluster indices
  sapply(1:nrow(nn), function(x) {
    # reorder the neighbors by size
    x <- nn[x, nn_order[x,]]
    
    x <- x[x <= nrow(data)]
    x <- object$cluster[x][x > 0][1]
    x[is.na(x)] <- 0L
    x
  })
  
}