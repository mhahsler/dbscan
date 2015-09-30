## This is the timing benchmark used by ELKI

#args <- commandArgs(trailingOnly = TRUE)
args <- list(file = "aloi-hsb-2x2x2.csv.gz", ncol = 8, minPts = 20, eps = 0.01)
print(as.data.frame(args))


library("microbenchmark")
library("dbscan")
requireNamespace("fpc")

rbind(
  format(packageDescription("dbscan", fields = c("Package", "Version"))),
  format(packageDescription("fpc", fields = c("Package", "Version")))
)

data <- read.csv(gzfile(as.character(args[1])), header = FALSE, sep = " ")
kdata <- data[,seq(1,as.integer(args[2]))]
dim(kdata)

#kdata <- kdata[1:10000,]

## compare results
d_dbscan <- dbscan::dbscan(kdata,
  eps = as.double(args[4]), minPts = as.integer(args[3]))

#d_fpc <- fpc::dbscan(kdata,
#  eps = as.double(args[4]), MinPts = as.integer(args[3]),
#  method="hybrid", seeds=FALSE)

load("d_fpc.rda")

rbind(
  dbscan = table(d_dbscan$cluster),
  fpc = table(d_fpc$cluster)
)

all.equal(d_dbscan$cluster, d_fpc$cluster)

## timing

t <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3])),
  times = 10
)

#system.time(
#  d_fpc <- fpc::dbscan(kdata, eps = as.double(args[4]), MinPts = as.integer(args[3]),
#    method="raw", seeds=FALSE)
#)
##  user   system  elapsed
## 2250.444    0.564 2248.322
#save(d_fpc, file="d_fpc.rda")

t_linear <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), search = "linear"),
  times = 1
)

t_bs50 <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), bucketSize = 50),
  times = 10
)

t_bs20 <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), bucketSize = 20),
  times = 10
)

t_bs5 <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), bucketSize = 5),
  times = 10
)

t_sr_std <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), splitRule = "STD"),
  times = 10
)

## dies
#t_sr_fair <- microbenchmark(
#  dbscan::dbscan(kdata,
#    eps = as.double(args[4]), minPts = as.integer(args[3]), splitRule = "FAIR"),
#  times = 10
#)

t_sr_midpt <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), splitRule = "MIDPT"),
  times = 10
)

t_approx <- microbenchmark(
  dbscan::dbscan(kdata,
    eps = as.double(args[4]), minPts = as.integer(args[3]), approx = .05),
  times = 10
)


time <- rbind(
  default = t,
  linear_search = t_linear,
  bucket_size_5 = t_bs5,
  bucket_size_20 = t_bs20,
  bucket_size_50 = t_bs50,
  splitting_rule_kd = t_sr_std,
  #splitting_rule_fair = t_sr_fair,
  splitting_rule_midpoint = t_sr_midpt,
  approx = t_approx
)

options(digits=3)

time$expr <- c(
  rep("default", 10),
  "linear_search",
  rep("bucket_size_5", 10),
  rep("bucket_size_20", 10),
  rep("bucket_size_50", 10),
  rep("splitting_rule_kd", 10),
  rep("splitting_rule_midpoint", 10),
  rep("approx", 10)
)

time

#Unit: seconds
#expr                      min     lq   mean median     uq    max neval
#approx                   1.42   1.45   1.49   1.45   1.47   1.74    10
#bucket_size_20           1.46   1.48   1.57   1.50   1.51   1.95    10
#bucket_size_5            1.68   1.69   1.80   1.74   1.92   2.19    10
#bucket_size_50           1.56   1.57   1.61   1.58   1.62   1.85    10
#default                  1.51   1.51   1.53   1.52   1.55   1.59    10
#linear_search          103.57 103.57 103.57 103.57 103.57 103.57     1
#splitting_rule_kd        1.54   1.54   1.56   1.55   1.59   1.63    10
#splitting_rule_midpoint  1.52   1.52   1.52   1.52   1.52   1.54    10

# fpc::dbscan takes 2250 seconds


## From ELKI maintainer (erich@vitavonni.de):
## Good results with indexes are 6-9 seconds in Java
## and 4 seconds in C code.
##
## Non-indexed good results are 220 seconds in Java
## Recent Weka is around 1169 seconds
## Octave around 950.
## Sklearn without index used to be around 1500 seconds, but now runs out of memory
## fpc version took 2931.62 (48min) seconds here
## dbscan package 6232 seconds (version 0.9-0).

