test_that("dbcv", {
  # From: https://github.com/FelSiq/DBCV
  #
  # Dataset	      MATLAB
  # dataset_1.txt	0.8576
  # dataset_2.txt	0.8103
  # dataset_3.txt	0.6319
  # dataset_4.txt	0.8688
  #
  # Original MATLAB implementation is at:
  #     https://github.com/pajaskowiak/dbcv/tree/main/data

  data(Dataset_1)
  x <- Dataset_1[, c("x", "y")]
  class <- Dataset_1$class
  #clplot(x, class)
  (db <- dbcv(x, class, metric = "sqeuclidean"))
  expect_equal(round(db$score, 2), 0.86)

  # detailed results from the Python implementation
  #dsc [0.00457826 0.00457826 0.0183068  0.0183068 ]
  #dspc [0.85627898 0.85627898 0.85627898 0.85627898]
  #vcs [0.99465331 0.99465331 0.97862052 0.97862052]
  #0.8575741400490697

  data(Dataset_2)
  x <- Dataset_2[, c("x", "y")]
  class <- Dataset_2$class
  #clplot(x, class)
  (db <- dbcv(x, class, metric = "sqeuclidean"))
  expect_equal(round(db$score, 2), 0.81)

  #dsc [19.06151967 15.6082 83.71522964 68.969]
  #dspc [860.2538 501.4376 501.4376 860.2538]
  #vcs [0.97784198 0.9688731  0.83304956 0.91982715]
  #0.8103343589093096

  # more data sets

  # data(Dataset_3)
  # x <- Dataset_3[, c("x", "y")]
  # class <- Dataset_3$class
  # #clplot(x, class)
  # (db <- dbcv(x, class, metric = "sqeuclidean"))
  #
  # data(Dataset_4)
  # x <- Dataset_4[, c("x", "y")]
  # class <- Dataset_4$class
  # #clplot(x, class)
  # (db <- dbcv(x, class, metric = "sqeuclidean"))

})
