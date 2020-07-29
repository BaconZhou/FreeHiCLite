library(FreeHiCLite)
library(hicrep)
create_matrix <- function(mat) {
    x <- as.integer(mat[,1] / 40000) + 1
    y <- as.integer(mat[,2] / 40000) + 1
    value <- mat[,3]
    require(Matrix)
    R <- max(max(x), max(y))
    m <- sparseMatrix(i = x, j = y, x = value, dims = c(R, R), symmetric = TRUE)
    as.matrix(m)
}


N <- 2000
maxBinX <- maxBinY <- 2000000L
maxCounts <- 10
binX <- sample(1:maxBinX, N, replace=TRUE)
binY <- sample(1:maxBinY, N, replace=TRUE)
counts <- sample(1:10, N, replace=TRUE)

seqDepth <- 1
countScale <- 0.5
noiseRate <- 0.5
neighborZeroRate <- 0
resolution <- 5000

## Matrix version
## matrix layout as x, y, counts
contacts <- matrix(0, N, 3)
contacts[,1] <- binX
contacts[,2] <- binY
contacts[,3] <- counts

contacts <- contacts[binY > binX, ]

sum(counts)

set.seed(123)
res1 <- FreeHiC(contacts, seqDepth, countScale, noiseRate = 0.1, neighborZeroRate, resolution)
set.seed(123)
res2 <- FreeHiC(contacts, seqDepth, countScale, noiseRate = 0.5, neighborZeroRate, resolution)

mc <- create_matrix(contacts)
m1 <- create_matrix(res1)
m2 <- create_matrix(res1)

get.scc(m1, mc, 40000, 1)$scc
get.scc(m1, m2, 40000, 1)$scc

head(res1)
head(res2)
dim(res1)
dim(res2)
