hullplot <- function(x, cl, col = NULL,
  cex = 0.5, hull_lwd = 1, hull_lty = 1, main = "Convex Cluster Hulls",
    solid=FALSE, alpha = .2, ...) {

  ### extract clustering (keep hierarchical OPTICSXi structure)
  if(is(cl, "optics") && !is.null(cl$clusters_xi)) clusters_xi <- cl
  else clusters_xi <- NULL

  if(is.list(cl)) cl <- cl$cluster
  if(!is.numeric(cl)) stop("could not get cluster assignment vector from cl.")

  #if(is.null(col)) col <- c("#000000FF", rainbow(n=max(cl)))
  if(is.null(col)) col <- palette()
  if(max(cl)+1L > length(col)) warning("not enought colors. some colors will be reused.")

  plot(x[,1:2], col = col[cl %% length(col) +1L], cex = cex, main = main, ...)

  col_poly <- adjustcolor(col, alpha.f = alpha)
  if(is.null(hull_lwd) || is.na(hull_lwd) || hull_lwd == 0) {
    hull_lwd <- 1
    border <- NA
  }else border <- NULL

  for(i in 1:max(cl)) {

    ### use all the points for OPTICSXi's hierarchical structure
    if(is.null(clusters_xi)) d <- x[cl==i,]
    else d <- x[clusters_xi$order[clusters_xi$clusters_xi$start[i] : clusters_xi$clusters_xi$end[i]],]

    ch <- chull(d)
    ch <- c(ch, ch[1])
    if(!solid)
      lines(d[ch,], col = col[i %% length(col) +1L], lwd=hull_lwd, lty=hull_lty)
    else
      polygon(d[ch,], col = col_poly[i %% length(col_poly) +1L],
        lwd=hull_lwd, lty=hull_lty, border = border)
  }
}
