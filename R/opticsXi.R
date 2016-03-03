#' WAMI.net - opticsXi function 
#' 
#' opticsXi is a hierarchical clustering algorithm implementation that was created as 
#' an extension of the OPTICS ordering algorithm. This is simple an implementation of the
#' clustering algorithm shown in Figure 19 of the OPTICS: Ordering Points To Identify the Clustering Structure 
#' paper published by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander. 
#' 
#' see http://fogo.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf for more details 
#' @param optics_scan optics object produced by an OPTICS ordering impelementation, such as 
#' the implementation from the package 'dbscan'. Required.   
#' @param xi The steepness value used to detect changes in relative cluster density 
#' @export
#' @examples
#' @export
opticsXi <- function(optics_scan, xi = 0.001)
{
  optics_scan$ord_rd <- optics_scan$reachdist[optics_scan$order]
  optics_scan$ixi <- (1 - xi)
  SetOfSteepDownAreas <- list()
  SetOfClusters <- list()
  index <- 1
  mib <- 0
  while (index < length(optics_scan$order))
  {
    mib <- ifelse(optics_scan$ord_rd[index] == Inf, mib, max(mib, optics_scan$ord_rd[index]))
    if (index == length(optics_scan$order))
      break
    if (DownPoint(index, optics_scan))
    {
   
      SetOfSteepDownAreas <- filterSetOfSteepDownAreas(index, SetOfSteepDownAreas, mib, optics_scan)
      D <- getSteepAreaBasic(i, "D", optics_scan)
      SetOfSteepDownAreas <- append(SetOfSteepDownAreas, list(D))
      index <- D$e + 1
      mib <- 0
      next 
    } else if (UpPoint(index,optics_scan))
    {
      U <- getSteepArea(i, "U", optics_scan)
      SetOfSteepDownAreas <- filterSetOfSteepDownAreas(index, SetOfSteepDownAreas, mib, optics_scan)
      index <- U$e + 1
      mib <- optics_scan$ord_rd[index]
      for (D in SetOfSteepDownAreas)
      {
        if (isValid(D, U, optics_scan) && satisfiesClusterConditions(D, U, optics_scan))
        {
          SetOfClusters <- append(SetOfClusters, list(computeClusterSE(D, U, optics_scan)))
        }
      }
    } else
    {
      index <- index + 1
    }
  }
  
  # Set up cluster variables
  clusterID <- 1; extClusterID <- 1
  optics_scan$clustersXi <- rep(0, length(optics_scan$order))
  cluster_df <- t(as.data.frame(SetOfClusters))
  
  # Get differences in point indices, effectively getting the size (number of points) in each cluster
  cluster_df <- cbind(cluster_df, cluster_df[, 2] - cluster_df[, 1])
  colnames(cluster_df) <- c("low", "high", "count")
  
  # Order clusters 
  cluster_df <- cluster_df[order(-cluster_df[, 3]), ]
  cluster_df_asc <- cluster_df[order(cluster_df[, 3]), ]
  optics_scan$clustersXiExtended <- rep(0, length(optics_scan$order))
  
  # Add clusters in order from biggest to smallest
  for (row in 1:nrow(cluster_df))
  {
    optics_scan$clustersXi[makeRange(cluster_df[row, ])] <- clusterID
    if (all(optics_scan$clustersXiExtended[makeRange(cluster_df_asc[row, ])] == 0))
    {
      optics_scan$clustersXiExtended[makeRange(cluster_df_asc[row, ])] <- extClusterID
      extClusterID <- extClusterID + 1
    }
    clusterID <- clusterID + 1
  }
  
  # Set class and return
  class(optics_scan) <- append(class(optics_scan), "opticsXi")
  return(optics_scan)
}
#' WAMI.net - plot.opticsXi function 
#'
#' opticsXi is a hierarchical clustering algorithm implementation that was created as 
#' an extension of the OPTICS ordering algorithm. This is simple an implementation of the
#' clustering algorithm shown in Figure 19 of the OPTICS: Ordering Points To Identify the Clustering Structure 
#' paper published by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander. 
#' 
#' see http://fogo.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf for more details 
#' @param optics_scan optics object produced by an OPTICS ordering impelementation, such as 
#' the implementation from the package 'dbscan'. Required.   
#' @param xi The steepness value used to detect changes in relative cluster density 
#' @export
#' @examples
#' @export
plot.opticsXi <- function(optics_results, min_clusters = F, point_data = NULL, points = F, polygons = F)
{
  colors <- rainbow(length(unique(optics_results$clustersXiExtended)))
  if (min_clusters == F && points == F)
  {
    plot(optics_results$ord_rd, type = "h", 
         col = optics_results$clustersXi,
         ylab = "Reachability dist.", xlab = "OPTICS order", 
         main = "Reachability Plot")
  } else if (min_clusters == T && points == F)
  {
    plot(optics_results$ord_rd, type = "h", 
         col = optics_results$clustersXiExtended, ylab = "Reachability dist.", 
         xlab = "OPTICS order", main = "Reachability Plot")
  } else if (min_clusters == T && points == T)
  {
    # Initial data plot we're looking at, black
    plot(point_data, pch = 20, cex = 0.4, col="black")
    
    # Add clusters incrementally as layers with random colors, ensuring the largest clusters get added first and the
    # smallest clusters last for maximum visibility
    min_cluster_ids <- sort(unique(optics_results$clustersXiExtended))[-1]
    for (clusterID in min_cluster_ids)
    {
      print(points(point_data[optics_results$order[which(optics_results$clustersXiExtended == clusterID)], ], col = sample(colours(), 1), pch = 20, cex = 0.85))
      print(paste("Plotting cluster:", clusterID, "of", length(min_cluster_ids)))
    }
  } else if (polygons == T)
  {
    #     for (location in location_clusters)
    #     {
    #       clocation <- location[chull(location$long, location$lat),]
    #       polygon(clocation$long, clocation$lat, border = "black", col = "green")
    #     }
  }
}




# Return T/F if the interval given by {D, U} satisfies cluster conditions 1, 2, 3a
satisfiesClusterConditions <- function(D, U, optics_scan)
{
  ## Cluster Conditions (ξ-cluster) Condition 1: DownAreaξ(D) ∧ s ∈ D Condition 2: UpAreaξ(U) ∧ e ∈ U Condition 3a: e – s ≥
  ## MinPts Note: In this implementation, if D is not an empty list, then DownAreaξ(D) = True (Similarly for UpAreaξ(U))
  return(length(D) > 0 && length(U) > 0 && U$e - D$s >= optics_scan$minPts)
}

computeClusterSE <- function(D, U, optics_scan)
{
  # Condition [a] (Default) (s, e) = (s_D, e_U) otherwise
  s <- D$s
  e <- U$e
  
  if (e == length(optics_scan$order))
  {
    return(c(s, e))
    # Condition [b] (s, e) = (max{ x ∈ D | r(x) > r(e_U + 1) }, e_U) if r(s_D) × (1 – ξ) ≥ r(e_U + 1)
  } else if (optics_scan$ord_rd[D$s] * optics_scan$ixi >= optics_scan$ord_rd[U$e + 1])
  {
    s <- max(Filter(function(x) optics_scan$ord_rd[x] > optics_scan$ord_rd[U$e + 1], D$s:D$e))
    e <- U$e
    # Condition [c] (s, e) = (s_D, min{ x ∈ U | r(x) < r(s_D) }) if r(e_U + 1) × (1 – ξ) ≥ r(s_D) [c]
  } else if (optics_scan$ord_rd[U$e + 1] * optics_scan$ixi >= optics_scan$ord_rd[D$s])
  {
    s <- D$s
    e <- min(Filter(function(x) optics_scan$ord_rd[x] < optics_scan$ord_rd[D$s], U$s:U$e))
  }
  return(c(s, e))
}

# We compare the end of the steep up area U multiplied by (1-ξ) to the mib-value of the steep down area D, thus
# satisfying (sc2*).  Ref (sc2*): max{ x | s_D < x < e_U} ≤ r(e_U) × (1 - ξ)
isValid <- function(D, U, optics_scan)
{
  # max(r[optics_order[(D$s+1):(U$e-1)]]) <= min(r[optics_order[c(D$s, U$e)]])
  return(D$mib <= optics_scan$ord_rd[U$e] * optics_scan$ixi)
}

# Whenever we encounter a new steep (up or down) area, we filter all steep down areas from SDASet whose start multiplied
# by (1-ξ) is smaller than the global mib-value, thus reducing the number of potential clusters significantly and
# satisfying (sc1*) at the same time (c.f. line (*) in figure 19).  
# Reference: (sc1*): max{ x | s_D < x < e_U ≤ r(s_D) × (1-ξ) }

filterSetOfSteepDownAreas <- function(i, SDASet, mib, optics_scan)
{
  SDASet <- lapply(SDASet, function(sda)
  {
    sda$mib <- max(optics_scan$ord_rd[(sda$e):i])
    sda
  })
  SDASet <- Filter(function(D) optics_scan$ord_rd[D$s] * optics_scan$ixi >= mib, SDASet)
  return(SDASet)
}


getSteepAreaBasic <- function (i, type, optics_scan)
{
  SteepPoint <- switch(type, D = DownPoint2, U = UpPoint2)
  NearlyFlat <- switch(type, D = nearlyDownPoint, U = nearlyUpPoint)
  I <- list(); s <- i; e <- i; scan_i <- 1

  while (e < length(optics_scan$order)-1)
  {
    if (SteepPoint(i))
    { 
      e <- e + 1 
      scan_i <- scan_i + 1
      next
    }
    if (!SteepPoint(i, optics_scan, 1.0) || scan_i - e > optics_scan$minPts)
    { break }
  }
  I$s <- s
  I$e <- e
  I$max <- optics_scan$ord_rd[s]
  I$mib <- 0
  return(I)
}
# An interval is called a ξ-steep upward area UpAreaξ(I) iff it satisfies the following conditions: s is a ξ-steep
# upward point: UpPointξ(s) e is a ξ-steep upward point: UpPointξ(e) ∀x, s < x ≤ e : r(x) ≥ r(x-1) each point between s
# and e is at least as high as its predecessor: I does not contain more than MinPts consecutive points that are not
# ξ-steep upward A ξ-steep downward area is defined analogously
#' @export
getSteepArea <- function(i, type, optics_scan)
{
  SteepPoint <- switch(type, D = DownPoint, U = UpPoint)
  NearlyFlat <- switch(type, D = nearlyDownPoint, U = nearlyUpPoint)
  I <- list()
  
  if (SteepPoint(i, optics_scan))
  {
    # Get the endpoints for the furthest steep area possible Ensure each point between s and e is at least as high as its
    # predecessor, also ensure the interval does not contain more than minPts consecutive non-steep-downward points
    s <- i
    e <- i
    cflat <- 0
    # Get the maximal area to be considered as a steep area 
    while (cflat < optics_scan$minPts && (SteepPoint(e, optics_scan) || NearlyFlat(e, optics_scan)))
    {
      e <- e + 1
      cflat <- ifelse(NearlyFlat(e, optics_scan), cflat + 1, 0)
      if (e == length(optics_scan$order)) 
        break
    }
    # Since the endpoints to a steep area are required to be steep points themselves, back up the right 
    # interval index until the last valid steep point was found
    while(SteepPoint(e, optics_scan) == F && e != length(optics_scan$order))
    {
      e <- e - 1 
    }
    I$s <- s
    I$e <- e
    I$mib <- 0
  }
  return(I)
}

# Wrapper for SteepUpPoint: Detect if the current index is the start of a steep up area
startOfSteepUpArea <- function(i, optics_scan)
{
  return(UpPoint(i,  optics_scan))
}

# Detect if the current index is the start of a steep down area
startOfSteepDownArea <- function(i, optics_scan)
{
  return(getSteepArea(i, "D", optics_scan))
}

# UpPointξ(p) p ⇔ r(p) ≤ r(p+1) × (1 – ξ)
UpPoint2 <- function(i, optics_scan, ixi)
{
  return(optics_scan$ord_rd[i] <= ixi * optics_scan$ord_rd[i + 1])
}

# Check if next point is still higher than predecessor, without steep ratio
nearlyUpPoint <- function(i, optics_scan)
{
  return(optics_scan$ord_rd[i] <= optics_scan$ord_rd[i + 1])
}

# DownPointξ(p) p ⇔ r(p) × (1 – ξ) ≤ r(p+1)
DownPoint2 <- function(i, optics_scan, ixi)
{
  return(optics_scan$ord_rd[i] * ixi >= optics_scan$ord_rd[i + 1])
}

# Check if next point is still lower than predecessor, without steep ratio
nearlyDownPoint <- function(i, optics_scan)
{
  return(optics_scan$ord_rd[i] >= optics_scan$ord_rd[i + 1])
}

makeRange <- function(v)
{
  return(seq(v[1], v[2]))
}

# Convert to a data frame cluster_df <- t(as.data.frame(clusters)) # Get differences in point indices, effectively
# getting the size (number of points) in each cluster cluster_df <- cbind(cluster_df, cluster_df[,2]-cluster_df[,1])
# colnames(cluster_df) <- c('low', 'high', 'count') # Order clusters by size (descending) cluster_df <-
# cluster_df[order(-cluster_df[,3]),]
getPointsInClusters <- function(clusters, centroids, optics_results)
{
  clusters_df <- t(as.data.frame(clusters))
  locs <- list()
  for (row in 1:nrow(clusters_df))
  {
    locs[[row]] <- as.data.frame(centroids[optics_results$order[clusters_df[row, ][1]:clusters_df[row, ][2]], ])
  }
  return(locs)
} 