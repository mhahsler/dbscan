library(dbscan)


x <- read.table("Work/data_DBCV/dataset_1.txt")
colnames(x) <- c("x", "y", "class")

cl <- x[, 3]
cl[cl < 0] <- 0
x[, 3] <- cl

plot(x[, 1:2], col = x[, 3] + 1L, asp = 1)

Dataset_1 <- x
save(Dataset_1, file="data/Dataset_1.rda", version = 2)

x <- read.table("Work/data_DBCV/dataset_2.txt")
colnames(x) <- c("x", "y", "class")

cl <- x[, 3]
cl[cl < 0] <- 0
x[, 3] <- cl

clplot(x[, 1:2], x[, 3])

Dataset_2 <- x
save(Dataset_2, file="data/Dataset_2.rda", version = 2)


x <- read.table("Work/data_DBCV/dataset_3.txt")
colnames(x) <- c("x", "y", "class")

cl <- x[, 3]
cl[cl < 0] <- 0
x[, 3] <- cl

clplot(x[, 1:2], x[, 3])

Dataset_3 <- x
save(Dataset_3, file="data/Dataset_3.rda", version = 2)

x <- read.table("Work/data_DBCV/dataset_4.txt")
colnames(x) <- c("x", "y", "class")

cl <- x[, 3]
cl[cl < 0] <- 0
x[, 3] <- cl

clplot(x[, 1:2], x[, 3])

Dataset_4 <- x
save(Dataset_4, file="data/Dataset_4.rda", version = 2)
