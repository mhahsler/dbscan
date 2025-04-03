# From: https://github.com/FelSiq/DBCV
#
# Dataset	Python (Scipy's Kruskal's)	Python (Translated MST algorithm)	MATLAB
# dataset_1.txt	0.8566	0.8576	0.8576
# dataset_2.txt	0.5405	0.8103	0.8103
# dataset_3.txt	0.6308	0.6319	0.6319
# dataset_4.txt	0.8456	0.8688	0.8688
#
# Original MATLAB implementation is at:
#     https://github.com/pajaskowiak/dbcv/tree/main/data


res <- c()

data(Dataset_1)
x <- Dataset_1[, c("x", "y")]
class <- Dataset_1$class
#clplot(x, class)
(db <- dbcv(x, class, metric = "sqeuclidean"))
res["ds1"] <- db$score


#dsc [0.00457826 0.00457826 0.0183068  0.0183068 ]
#dspc [0.85627898 0.85627898 0.85627898 0.85627898]
#vcs [0.99465331 0.99465331 0.97862052 0.97862052]
#0.8575741400490697

data(Dataset_2)
x <- Dataset_2[, c("x", "y")]
class <- Dataset_2$class
#clplot(x, class)
(db <- dbcv(x, class, metric = "sqeuclidean"))
res["ds2"] <- db$score

#dsc [19.06151967 15.6082     83.71522964 68.969     ]
#dspc [860.2538 501.4376 501.4376 860.2538]
#vcs [0.97784198 0.9688731  0.83304956 0.91982715]
#0.8103343589093096


data(Dataset_3)
x <- Dataset_3[, c("x", "y")]
class <- Dataset_3$class
#clplot(x, class)
(db <- dbcv(x, class, metric = "sqeuclidean"))
res["ds3"] <- db$score

data(Dataset_4)
x <- Dataset_4[, c("x", "y")]
class <- Dataset_4$class
#clplot(x, class)
(db <- dbcv(x, class, metric = "sqeuclidean"))
res["ds4"] <- db$score

cbind(dbscan = round(res, 2), MATLAB = c(0.85, 0.81, 0.63, 0.87))

