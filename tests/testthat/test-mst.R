test_that("mst", {
  draw_mst <- function(x, m) {
    plot(x)
    text(x, labels = 1:nrow(x), pos = 1)
    for (i in seq(nrow(m))) {
      from_to <- rbind(x[m[i, 1], ], x[m[i, 2], ])
      lines(from_to[, 1], from_to[, 2])
    }
  }

  x <- rbind(c(0, 0), c(0, 1), c(1, 1))
  d <- dist(x)
  (m <- mst(d, n = nrow(x)))

  draw_mst(x, m)

  expect_equal(m, structure(
    c(2, 3, 1, 2, 1, 1),
    dim = 2:3,
    dimnames = list(NULL, c("from", "to", "weight"))
  ))

  x <- rbind(c(0, 0),
             c(1, 0),
             c(0, 1),
             c(1, 1),
             c(2, 1),
             c(1, 2),
             c(.7, 1),
             c(.7, .7),
             c(.7, 1.3))
  d <- dist(x)
  (m <- mst(d, n = nrow(x)))

  draw_mst(x, m)

  expect_equal(m, structure(
    c(
      2,
      3,
      4,
      5,
      6,
      7,
      8,
      9,
      8,
      7,
      7,
      4,
      9,
      8,
      1,
      7,
      0.761577310586391,
      0.7,
      0.3,
      1,
      0.761577310586391,
      0.3,
      0.989949493661166,
      0.3
    ),
    dim = c(8L, 3L),
    dimnames = list(NULL, c("from", "to", "weight"))
  ))
})

test_that("dist_subset", {
  x <- rbind(c(0, 0),
             c(1, 0),
             c(0, 1),
             c(1, 1),
             c(2, 1),
             c(1, 2),
             c(.7, 1),
             c(.7, .7),
             c(.7, 1.3))
  d <- dist(x)
  m <- as.matrix(d)

  s <- c(1:3, 6)
  (d_sub <- dist_subset(d, s))
  (m_sub <- m[s,s])

  expect_equal(unname(as.matrix(d_sub)), unname(m_sub))
})
