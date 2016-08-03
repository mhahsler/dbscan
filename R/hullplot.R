hullplot <- function(x, cl, col = NULL,
  cex = 0.5, hull_lwd = 1, hull_lty = 1, main = "Convex Cluster Hulls",
    solid=FALSE, alpha = .2, ...) {

  if(is.list(cl)) cl <- cl$cluster
  if(!is.numeric(cl)) stop("could not get cluster assignment vector from cl.")

  if(is.null(col)) col <- palette()
  if(max(cl)+1L > length(col)) warning("not enought colors. some colors will be reused.")


  plot(x[,1:2], col = col[cl %% length(col) +1L], cex = cex, main = main, ...)

  col_poly <- adjustcolor(col, alpha.f = alpha)
  if(is.null(hull_lwd) || is.na(hull_lwd) || hull_lwd == 0) {
    hull_lwd <- 1
    border <- NA
  }else border <- NULL

  for(i in 1:max(cl)) {
    d <- x[cl==i,]
    ch <- chull(d)
    ch <- c(ch, ch[1])
    if(!solid)
      lines(d[ch,], col = col[i %% length(col) +1L], lwd=hull_lwd, lty=hull_lty)
    else
      polygon(d[ch,], col = col_poly[i %% length(col_poly) +1L],
        lwd=hull_lwd, lty=hull_lty, border = border)
  }
}
