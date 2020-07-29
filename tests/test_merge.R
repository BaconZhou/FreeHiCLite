library(dplyr)
library(purrr)
library(tidyr)

N = 1000

x <- random_matrix(N)
y <- random_matrix(N)
z <- random_matrix(N)
w <- random_matrix(N)
q <- random_matrix(N)

tmp <- mergeContactMatrixs(x, y, z, w, q)
tmp2 <- summaryContactMatrixs(tmp)
summary(tmp2[,3])
