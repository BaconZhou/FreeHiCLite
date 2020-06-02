#' @title 
#' 
#' FreeHiC simulator
#' 
#' @description 
#' 
#' This function performs FreeHiC based on contact matrix only.
#' 
#' @param contacts A \code{\link[base]{list}} of contacts matrix. The list element should named by chromosome pair.
#' \itemize{
#' \item{\code{\link[base]{matrix}}}{A matrix with columns (x, y, counts).}
#' \item{\code{\link[base]{data.frame}}}{A data.frame with columns (chr1, x, y, chr2, counts).}
#' \item{\code{\link[base]{list}}}{A list of matrix, where the matrix has columns (x, y, counts).}
#' }
#' where \code{x} stands for the first chromosome location, \code{y} stands for the second chromosome location and \code{counts} is the interaction counts.
#' @param seqDepth The desired sequence depth.
#' @param countScale The scale of counts. A number larger than 0. If both \code{seqDepth} and \code{countScale} are provided. Choose the larger one.
#' @param noiseRate The noise rate for contact matrix. 0 - 1 scale
#' @param neighborZeroRate The rate for neighborhood noise rate. 0 - 1 scale
#' @param resolution The resolution used in the contacts matrix. A positive number.
#' 
#' @return A list or matrix or data.frame with the same format as \code{contacts}
#' 
#' @examples 
#' 
#' library(FreeHiCLite)
#' 
#' N <- 2000
#' maxBinX <- maxBinY <- 2000000L
#' maxCounts <- 10
#' binX <- sample(1:maxBinX, N, replace=TRUE)
#' binY <- sample(1:maxBinY, N, replace=TRUE)
#' counts <- sample(1:10, N, replace=TRUE)
#' 
#' seqDepth <- 20000
#' countScale <- 0
#' noiseRate <- 0.1
#' neighborZeroRate <- 0
#' resolution <- 5000
#' 
#' ## Matrix version
#' ## matrix layout as x, y, counts
#' contacts <- matrix(0, N, 3)
#' contacts[,1] <- binX
#' contacts[,2] <- binY
#' contacts[,3] <- counts
#' 
#' res <- FreeHiC(contacts, seqDepth, countScale, noiseRate, neighborZeroRate, resolution)
#' head(res)
#' 
#' ## List version
#' N2 = 3000
#' binX2 <- sample(1:maxBinX, N2, replace=TRUE)
#' binY2 <- sample(1:maxBinY, N2, replace=TRUE)
#' counts2 <- sample(1:10, N2, replace=TRUE)
#' 
#' contacts2 <- matrix(0, N2, 3)
#' contacts2[,1] <- binX2
#' contacts2[,2] <- binY2
#' contacts2[,3] <- counts2
#' 
#' contactsMap <- list("1_1" = contacts, "2_2" = contacts2)
#' res <- FreeHiC(contactsMap, seqDepth, countScale, noiseRate, neighborZeroRate, resolution)
#' str(res)
#' 
#' ## Dataframe version
#' 
#' chr1 <- c(rep("1", N), rep("2", N2))
#' chr2 <- c(rep("1", N), rep("2", N2))
#' 
#' contactsAll <- rbind(contacts, contacts2)
#' 
#' contactsDf <- data.frame(chr1 = chr1, x = contactsAll[,1], chr2 = chr2, y = contactsAll[,2], counts = contactsAll[,3])
#' 
#' res <- FreeHiC(contactsDf, seqDepth, countScale, noiseRate, neighborZeroRate, resolution)
#' head(res)
#' 
#' @references 
#' Zheng, Ye, Keles, Sunduz FreeHi-C simulates high-fidelity Hi-C data for benchmarking and data augmentation. 
#' \emph{Nature Methods} \strong{17}, 37–40 (2020). \href{https://doi.org/10.1038/s41592-019-0624-3}{doi}
#' 
#' @importFrom methods is
#' @export
FreeHiC <- function(contacts, seqDepth = NULL, countScale = 1,
                    noiseRate = 0.0, neighborZeroRate = 0.0,  resolution = 50000L) {
  
  stopifnot(is.null(seqDepth) | seqDepth >= 0)
  stopifnot(countScale > -1)
  stopifnot(noiseRate >= 0 & noiseRate < 1)
  stopifnot(neighborZeroRate >= 0 & neighborZeroRate < 1)
  stopifnot(resolution >= 0)
  
  if (is(contacts, "list")) {
    ans <- .FreeHiCList(contactsMap = contacts, seqDepth=seqDepth, countScale = countScale,
                         noiseRate = noiseRate, neighborZeroRate = neighborZeroRate, resolution=resolution)
  } else if (is(contacts, "data.frame")) {
    ans <- .FreeHiCDf(contactsDf = contacts, seqDepth=seqDepth, countScale = countScale,
                         noiseRate = noiseRate, neighborZeroRate = neighborZeroRate, resolution=resolution)
  } else if (is(contacts, "matrix")) {
    ans <- .FreeHiCMatrix(contactsMatrix = contacts, seqDepth=seqDepth, countScale = countScale,
                       noiseRate = noiseRate, neighborZeroRate = neighborZeroRate, resolution=resolution)
  } else {
    stop("Not implemented. Check ?FreeHIC")
  }
  return(ans)
}

.FreeHiCList <- function(contactsMap, seqDepth, countScale,
                         noiseRate, neighborZeroRate,resolution) {
  pairs <- names(contactsMap)
  stopifnot(!is.null(pairs))
  return(hicDataSimuList(contactRecords = contactsMap, names = pairs, resolution = resolution, 
                         sequenceDepth = seqDepth, countScale = countScale, 
                         noiseRate= noiseRate, neighborZeroRate = neighborZeroRate))
}

.FreeHiCDf <- function(contactsDf, seqDepth, countScale,
                           noiseRate, neighborZeroRate,resolution) {
  stopifnot(NCOL(contactsDf) == 5)
  contactsMap <- .dfToList(contactsDf)
  resList <- .FreeHiCList(contactsMap = contactsMap, 
                          seqDepth = seqDepth, countScale = countScale,
                          noiseRate = noiseRate, neighborZeroRate = neighborZeroRate, resolution = resolution)
  pairs <- names(resList)
  counts <- sapply(resList, NROW)
  chrs <- sapply(pairs, FUN = function(x) {strsplit(x, "_")})
  
  chr1 <- c()
  chr2 <- c()
  
  for (i in seq_along(counts)) {
    chr1 <- c(chr1, rep(chrs[[i]][1], counts[i]))
    chr2 <- c(chr2, rep(chrs[[i]][2], counts[i]))
  }
  
  mat <- do.call('rbind', resList)
  rm(resList)
  
  df <- data.frame(chr1 = chr1, x = mat[,1], chr2 = chr2, y = mat[,2], counts = mat[,3])
  names(df) <- names(contactsDf)
  return (df)
}
  
.FreeHiCMatrix <- function(contactsMatrix, seqDepth, countScale,
                       noiseRate, neighborZeroRate,resolution) {
  stopifnot(NCOL(contactsMatrix) == 3)
  return(hicDataSimuMatrix(contactRecords = contactsMatrix, 
                    resolution = resolution, 
                    sequenceDepth = seqDepth, countScale = countScale,
                    noiseRate = noiseRate, neighborZeroRate = neighborZeroRate))
}

#' @title
#' 
#' FreeHiC directly from hic file
#' 
#' @param file Filename can be a local path or remote path. 
#' The remote path full list can be obtained from \url{http://aidenlab.org/data.html}.
#' @param chromosomes A vector contains all the chromosomes. 
#' For example c("chr1", "chr2"), the resulting contact matrixs will include all the pairs of interaction (chr1_chr1, chr1_chr2, chr2_chr2). 
#' @param pairs A vector contains all the pair. 
#' The pair take format as "1_1" or "chr1_chr1", both means the contact between chromosome1 and chromosome1. If \code{pair} presents, \code{chromosomes} argument will be ignore.
#' @param unit Unit only supports c("BP", "FRAG"). 
#' "BP" means base-pair, and "FRAG" means fragment.
#' @param resolution The resolution used in the contacts matrix. A positive number.
#' The resolution must be a value from following list.
#' \itemize{
#' \item{\strong{unit: BP}}
#' \itemize{
#' \item 2500000
#' \item 1000000
#' \item 500000
#' \item 250000
#' \item 100000
#' \item 50000
#' \item 25000
#' \item 10000
#' \item 5000
#' }
#' \item{\strong{unit: FRAG}}
#' \itemize{
#' \item 500
#' \item 250
#' \item 100
#' \item 50
#' \item 20
#' \item 5
#' \item 2
#' \item 1
#' }
#' }
#' @param seqDepth The desired sequence depth.
#' @param countScale The scale of counts. A number larger than 0. If both \code{seqDepth} and \code{countScale} are provided. Choose the larger one.
#' @param noiseRate The noise rate for contact matrix. 0 - 1 scale
#' @param neighborZeroRate The rate for neighborhood noise rate. 0 - 1 scale
#' 
#' @return A list contains simulated contact matrix and basic information of hic data
#' 
#' @examples
#' 
#' library(FreeHiCLite)
#' 
#' ## Local file location
#' localFilePath = system.file("extdata", "example.hic", package = "FreeHiCLite")
#' 
#' ## Chromosomes needs to be extra
#' chromosomes = c("chr1", "chr2")
#' 
#' ## Pairs needs to be extra
#' pairs = c("1_1", "1_2")
#' unit = "BP"
#' resolution = 500000L
#' 
#' seqDepth <- 20000
#' countScale <- 0
#' noiseRate <- 0.1
#' neighborZeroRate <- 0
#' 
#' 
#' ## pass chrosomes into function, it will contains all the interaction pairs
#' 
#' res <- FreeHiCJuicer(file=localFilePath, chromosomes=chromosomes, pairs = NULL, unit=unit, resolution=resolution,
#' seqDepth = seqDepth, countScale = countScale, noiseRate = noiseRate, neighborZeroRate = neighborZeroRate)
#' 
#' str(res)
#' 
#' 
#' @export
FreeHiCJuicer <- function(file, chromosomes = NULL, pairs = NULL, 
                          unit = c("BP", "FRAG"), seqDepth = NULL, countScale = 1,
                          noiseRate = 0.0, neighborZeroRate = 0.0, resolution = 50000L  ) {
  defResolution <- list("BP" = c(2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000),
                        "FRAG" = c(500, 250, 100, 50, 20, 5, 2, 1))
  isHttp <- startsWith(file, 'http')
  if (!isHttp) {
    file = normalizePath(file)
    stopifnot(file.exists(file))
  } 
  message("You are trying to access file from: ", file, ". Make sure you provide a valid path.")
  
  if (isHttp) message("Remote file can take a while.\nThe full list of remote files can be found in http://aidenlab.org/data.html")
  
  ## process chromsome
  
  chromosomesChar <- NULL
  if (!is.null(chromosomes)) {
    chromosomesChar <- .chromosomes_clean(chromosomes)
  } 
  
  ## process pairs
  pairClean <- NULL
  if (!is.null(pairs)) {
    pairClean <- .parse_pair(pairs)
    if (is.logical(pairClean)) {stop("check input pair! ", paste(pairs, collapse = ", "))}
  } else {
    pairClean <- .chrosomes_to_pair(chromosomesChar)
  }
  message("Use pairs: ", paste(pairClean, collapse = ", "))
  
  unit = match.arg(unit)
  
  if (length(resolution) > 1) {
    stop("You provide resolution: ", paste0(resolution, collapse = ','), ". Please change the resolution to a number.\n",
         "From resolution (", unit, "): ", paste0(defResolution[[unit]], collapse = ','))
  }
  
  hicInfo <- readJuicerInformation(file)
  
  {
    unitResolution <- hicInfo[['resolution']]
    
    stopifnot(unit %in% names(unitResolution))
    stopifnot(resolution %in% unitResolution[[unit]])
    
    chrAll <- hicInfo[["chromosomeSizes"]][['chromosome']]
    .check_chr(pairClean, chrAll)
  }
  
  ans <- hicDataSimuFromFile(file, isHttp, pairClean, unit, resolution, 
                             seqDepth, countScale, noiseRate, neighborZeroRate)
  
  ans$information <- hicInfo
  
  return (ans)
}

.dfToList <- function(contactsDf) {
  colnames(contactsDf) = c("chr1", "x", "chr2", "y", "counts")
  pairs <- paste(contactsDf$chr1, contactsDf$chr2, sep="_")
  return(lapply(split(contactsDf[,c("x", "y", "counts")], pairs), as.matrix))
}

