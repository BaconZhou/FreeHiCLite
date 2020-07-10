mergeContactMatrixs <- function(..., normalized = TRUE) {
    dats <- list(...)
    N <- length(dats)
    if (N == 1) return(dats[[1]])
    for (i in 1:N) {
        colnames(dats[[i]]) <- c("x", "y", paste0("counts.", i))
    }
    mergeMatrix <- Reduce(function(d1, d2) merge(d1, d2, by = c("x", "y"), all = TRUE), dats)   
    rm(dats)
    
    if (normalized) {
        nc <- NCOL(mergeMatrix)
        for(j in 3:nc) mergeMatrix[is.na(mergeMatrix[,j]),j] = 0
        seqDep <- colSums(mergeMatrix[,3:nc])
        maxDep <- max(seqDep)
        for(j in 3:nc) mergeMatrix[,j] <- mergeMatrix[,j] / seqDep[j - 2] * maxDep 
    }
    return (mergeMatrix)
}

summaryContactMatrixs <- function(mergeMatrix) {
    nc <- NCOL(mergeMatrix)
    stopifnot(nc >= 3)
    countMedian <- apply(mergeMatrix[, 3:nc], 1, median, na.rm=TRUE)
    countMedian[countMedian < 1] = 0.45
    return(cbind(mergeMatrix[,1:2], countMedian))
}
