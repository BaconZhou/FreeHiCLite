#' @title
#' 
#' FreeHiC spikein
#' 
#' @description 
#' 
#' Add spikein to contact matrix.
#' 
#' @param contactBackground The contact matrix. It can be 
#' \itemize{
#' \item{\code{\link[base]{matrix}}}{A matrix with columns (x, y, counts).}
#' \item{\code{\link[base]{data.frame}}}{A data.frame with columns (chr1, x, y, chr2, counts).}
#' \item{\code{\link[base]{list}}}{A list of matrix, where the matrix has columns (x, y, counts).}
#' }
#' where \code{x} stands for the first chromosome location, \code{y} stands for the second chromosome location and \code{counts} is the interaction counts.
#' @param contactSpikeInSignal The spikein signal. It should has the exact same format as \code{contactBackground}
#' @param kernelSmooth TRUE or FALSE. Whether to perform kernel smoothing.
#' @param bandwidth The bandwidth used in the kernel smooth.
#' 
#' @return A same format as \code{contactBackground}
#' 
#' @details 
#' 
#' Spikein will add signals to the background. Also use Gaussian kernel smooth with bandwidth.
#' 
#' @examples 
#' 
#' library(FreeHiCLite)
#' 
#' N <- 2000
#' Ns <- 200
#' maxBinX <- maxBinY <- 2000000L
#' maxCounts <- 10
#' binX <- sample(1:maxBinX, N, replace=TRUE)
#' binY <- sample(1:maxBinY, N, replace=TRUE)
#' counts <- sample(1:10, N, replace=TRUE)
#' 
#' kernelSmooth = TRUE
#' bandwidth = 50000L
#' 
#' 
#' ## Matrix version
#' ## matrix layout as x, y, counts
#' 
#' contacts <- matrix(0, N, 3)
#' contacts[,1] <- binX
#' contacts[,2] <- binY
#' contacts[,3] <- counts
#' 
#' spikeIn <- contacts[sample(1:N, Ns),]
#' hist(spikeIn[,3])
#' 
#' spikeIn[,3] <- spikeIn[,3] * sample(seq(0, 10, 0.5), Ns, replace=TRUE)
#' hist(spikeIn[,3])
#' 
#' res <- FreeSpikeIn(contacts, spikeIn, kernelSmooth = kernelSmooth, bandwidth = bandwidth)
#' head(res)
#' 
#' ## List version
#' N2 = 3000
#' Ns2 = 300
#' binX2 <- sample(1:maxBinX, N2, replace=TRUE)
#' binY2 <- sample(1:maxBinY, N2, replace=TRUE)
#' counts2 <- sample(1:10, N2, replace=TRUE)
#' 
#' contacts2 <- matrix(0, N2, 3)
#' contacts2[,1] <- binX2
#' contacts2[,2] <- binY2
#' contacts2[,3] <- counts2
#' 
#' spikeIn2 <- contacts2[sample(1:N2, Ns2),]
#' spikeIn2[,3] <- spikeIn2[,3] * sample(seq(0, 10, 0.5), Ns2, replace=TRUE)
#' 
#' contactsBackgroup <- list('1_1' = contacts, '2_2' = contacts2)
#' spikeInlist <- list('1_1' = spikeIn, '2_2' = spikeIn2)
#' 
#' res <- FreeSpikeIn(contactsBackgroup, spikeInlist, 
#' kernelSmooth = kernelSmooth, bandwidth = bandwidth)
#' str(res)
#' 
#' ## Dataframe version
#' 
#' chr1 <- c(rep('1', N), rep('2', N2))
#' chr2 <- c(rep('1', N), rep('2', N2))
#' 
#' schr1 <- c(rep('1', Ns), rep('2', Ns2))
#' schr2 <- c(rep('1', Ns), rep('2', Ns2))
#' 
#' contactsAll <- rbind(contacts, contacts2)
#' spikeInAll <- rbind(spikeIn, spikeIn2)
#' 
#' contactsDf <- data.frame(chr1 = chr1, x = contactsAll[,1], 
#' chr2 = chr2, y = contactsAll[,2], counts = contactsAll[,3])
#' spikeInDf <- data.frame(chr1 = schr1, x = spikeInAll[,1], 
#' chr2 = schr2, y = spikeInAll[,2], counts = spikeInAll[,3])
#' res <- FreeSpikeIn(contactsDf, spikeInDf, 
#' kernelSmooth = kernelSmooth, bandwidth = bandwidth)
#' head(res)
#' 
#' 
#' @importFrom methods is
#' 
#' @export
FreeSpikeIn <- function(contactBackground, contactSpikeInSignal, kernelSmooth = TRUE, bandwidth = 500000L) {
    
    stopifnot(is.logical(kernelSmooth))
    stopifnot(bandwidth > 0)
    if (methods::is(contactBackground, "list")) {
        stopifnot(methods::is(contactSpikeInSignal, "list"))
        ans <- .FreeSpikeInList(contactBackground = contactBackground, contactSpikeInSignal = contactSpikeInSignal, 
            kernelSmooth = kernelSmooth, bandwidth = bandwidth)
    } else if (methods::is(contactBackground, "data.frame")) {
        stopifnot(methods::is(contactSpikeInSignal, "data.frame"))
        ans <- .FreeSpikeInDf(contactBackground = contactBackground, contactSpikeInSignal = contactSpikeInSignal, 
            kernelSmooth = kernelSmooth, bandwidth = bandwidth)
    } else if (methods::is(contactBackground, "matrix")) {
        stopifnot(methods::is(contactSpikeInSignal, "matrix"))
        ans <- .FreeSpikeInMatrix(contactBackground = contactBackground, contactSpikeInSignal = contactSpikeInSignal, 
            kernelSmooth = kernelSmooth, bandwidth = bandwidth)
    } else {
        stop("Not implemented. Check ?FreeSpikeIn")
    }
    return(ans)
}

.FreeSpikeInList <- function(contactBackground, contactSpikeInSignal, kernelSmooth, bandwidth) {
    pairs1 <- sort(names(contactBackground))
    pairs2 <- sort(names(contactSpikeInSignal))
    stopifnot(all(pairs1 == pairs2))
    
    ans <- list()
    for (i in seq_along(pairs1)) {
        pair <- pairs1[i]
        ans[[pair]] <- spikein(contactBackground[[pair]], contactSpikeInSignal[[pair]], bandwidth = bandwidth, 
            smooth = kernelSmooth)
    }
    return(ans)
}

.FreeSpikeInDf <- function(contactBackground, contactSpikeInSignal, kernelSmooth, bandwidth) {
    
    stopifnot(NCOL(contactBackground) == 5)
    contactsMap <- .dfToList(contactBackground)
    spikeInMap <- .dfToList(contactSpikeInSignal)
    resList <- .FreeSpikeInList(contactBackground = contactsMap, contactSpikeInSignal = spikeInMap, 
        kernelSmooth = kernelSmooth, bandwidth = bandwidth)
    pairs <- names(resList)
    counts <- sapply(resList, NROW)
    chrs <- sapply(pairs, FUN = function(x) {
        strsplit(x, "_")
    })
    
    chr1 <- c()
    chr2 <- c()
    
    for (i in seq_along(counts)) {
        chr1 <- c(chr1, rep(chrs[[i]][1], counts[i]))
        chr2 <- c(chr2, rep(chrs[[i]][2], counts[i]))
    }
    
    mat <- do.call("rbind", resList)
    rm(resList)
    
    df <- data.frame(chr1 = chr1, x = mat[, 1], chr2 = chr2, y = mat[, 2], counts = mat[, 3])
    names(df) <- names(contactBackground)
    return(df)
    
}
.FreeSpikeInMatrix <- function(contactBackground, contactSpikeInSignal, kernelSmooth, bandwidth) {
    
    stopifnot(NCOL(contactBackground) == 3)
    stopifnot(NCOL(contactSpikeInSignal) == 3)
    
    return(spikein(contactBackground, contactSpikeInSignal, bandwidth = bandwidth, smooth = kernelSmooth))
    
}
