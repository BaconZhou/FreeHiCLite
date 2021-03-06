#' @title 
#' Juicebox hic file reader
#' 
#' @description 
#' Read \href{https://github.com/aidenlab/juicer/wiki/Data}{.hic} file generated by \href{https://github.com/aidenlab/juicebox}{Juicebox}. 
#' Currently can be used both local and remote file.
#' 
#' @details 
#' This function is heavily adopted from both java version juicebox 
#' \href{https://github.com/aidenlab/juicer/wiki/Data-Extraction}{Dump} function 
#' and c++ version \href{https://github.com/aidenlab/straw}{straw}.
#' 
#' @param file Filename can be a local path or remote path. 
#' The remote path full list can be obtained from \url{http://aidenlab.org/data.html}.
#' @param chromosomes A vector contains all the chromosomes. 
#' For example c('chr1', 'chr2'), the resulting contact matrixes will include all the pairs of interaction (chr1_chr1, chr1_chr2, chr2_chr2). 
#' @param pairs A vector contains all the pair. 
#' The pair take format as '1_1' or 'chr1_chr1', both means the contact between chromosome1 and chromosome1. If \code{pairs} presents, \code{chromosomes} argument will be ignore.
#' @param unit Unit only supports c('BP', 'FRAG'). 
#' 'BP' means base-pair, and 'FRAG' means fragment.
#' @param resolution The desired resolution of the contact matrix. 
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
#' @param verbose TRUE or FALSE. Whether print information.
#' 
#' @return A list object includes following items.
#' \item{\code{contact}}{A list contains all the contact matrix. The keys of list is coded as \code{chrA_chrB}.}
#' \item{\code{information}}{A list contains basic information of the contact matrix.}
#' \itemize{
#' \item{\code{genomeID}}{The genome id of the current hic file.}
#' \item{\code{resolution}}{The list of current hic file available resolution.}
#' \item{\code{chromosomeSizes}}{A dataframe of chromosome size (chromosome, size). Can be used in juicer pre function for different genome.}
#' }
#' \item{\code{settings}}{A list contains file name, unit, resolution, and chromosomes.}
#' 
#' @examples
#' 
#' library(FreeHiCLite)
#' 
#' ## Remote file location. The reomte file include downloading, it may take a while
#' remoteFilePath = 'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic'
#' 
#' ## Local file location
#' localFilePath = system.file('extdata', 'example.hic', package = 'FreeHiCLite')
#' 
#' ## Chromosomes needs to be extracted
#' chromosomes = c('chr1', 'chr2')
#' 
#' ## Pairs needs to be extracted
#' pairs = c('1_1', '1_2')
#' unit = 'BP'
#' resolution = 500000L
#' 
#' ## pass chrosomes into function, it will contains all the interaction pairs
#' 
#' dat <- readJuicer(file=localFilePath, chromosomes=chromosomes, 
#' pairs = NULL, unit=unit, resolution=resolution)
#' 
#' print(names(dat[['contact']]))
#' 
#' ## pass pairs into function, it will contains only the given pairs
#' start = Sys.time()
#' dat <- readJuicer(file=localFilePath, chromosomes=NULL, 
#' pairs = pairs, unit=unit, resolution=resolution)
#' end = Sys.time()
#' 
#' print(end - start)
#' str(dat) 
#' 
#' ## Access each contact matrix
#' 
#' contacts = dat[['contact']]
#' 
#' ## chromosome 1 vs chromosome 1
#' key = '1_1'
#' contactMatrix <- contacts[[key]]
#' head(contactMatrix)
#' 
#' ## chromosome 1 vs chromosome 2
#' key = '1_2'
#' contactMatrix <- contacts[[key]]
#' head(contactMatrix)
#' 
#' 
#' 
#' @export
readJuicer <- function(file, chromosomes = NULL, pairs = NULL, unit = c("BP", "FRAG"), resolution, 
    verbose = FALSE) {
    defResolution <- list(BP = c(2500000, 1e+06, 5e+05, 250000, 1e+05, 50000, 25000, 10000, 5000), 
        FRAG = c(500, 250, 100, 50, 20, 5, 2, 1))
    isHttp <- startsWith(file, "http")
    if (!isHttp) {
        file = normalizePath(file)
        stopifnot(file.exists(file))
    }
    if (verbose) 
        message("You are trying to access file from: ", file, ". Make sure you provide a valid path.")
    
    if (isHttp & verbose) 
        message("Remote file can take a while.\nThe full list of remote files can be found in http://aidenlab.org/data.html")
    
    ## process chromosome
    
    chromosomesChar <- NULL
    if (!is.null(chromosomes)) {
        chromosomesChar <- .chromosomes_clean(chromosomes)
    }
    
    ## process pairs
    pairClean <- NULL
    if (!is.null(pairs)) {
        pairClean <- .parse_pair(pairs)
        if (is.logical(pairClean)) {
            stop("check input pair! ", paste(pairs, collapse = ", "))
        }
    } else {
        pairClean <- .chrosomes_to_pair(chromosomesChar)
    }
    if (verbose) 
        message("Use pairs: ", paste(pairClean, collapse = ", "))
    
    unit = match.arg(unit)
    
    if (length(resolution) > 1) {
        stop("You provide resolution: ", paste0(resolution, collapse = ","), ". Please change the resolution to a number.\n", 
            "From resolution (", unit, "): ", paste0(defResolution[[unit]], collapse = ","))
    }
    
    hicInfo <- readJuicerInformation(file, verbose = verbose)
    
    {
        unitResolution <- hicInfo[["resolution"]]
        
        stopifnot(unit %in% names(unitResolution))
        stopifnot(resolution %in% unitResolution[[unit]])
        
        chrAll <- hicInfo[["chromosomeSizes"]][["chromosome"]]
        .check_chr(pairClean, chrAll)
    }
    
    ans <- hicDataExtra(fileName = file, isHttp = isHttp, pair = pairClean, unit = unit, resolution = resolution)
    {
        cs <- ans$information$chromosomeSizes
        ord <- order(cs[["chromosome"]])
        cs <- cs[ord, ]
        ans$information$chromosomeSizes <- cs
    }
    
    ans$settings <- list(unit = unit, resolution = resolution, chromosomes = chromosomes, file = file)
    class(ans) = c("juicer", "freehic")
    ans
}



#' @title 
#' Juicebox hic file information
#' 
#' @description 
#' Read \href{https://github.com/aidenlab/juicer/wiki/Data}{.hic} file generated by \href{https://github.com/aidenlab/juicebox}{Juicebox}. Currently can be used for both local and remote file.
#' 
#' @details 
#' 
#' This function is heavily adopted from both java version juicebox \href{https://github.com/aidenlab/juicer/wiki/Data-Extraction}{Dump} function 
#' and c++ version \href{https://github.com/aidenlab/straw}{straw}.
#' 
#' In case people needs to use juice box \href{https://github.com/aidenlab/juicer/wiki/Pre}{Pre} to generate hic file. 
#' Check \code{\link[FreeHiCLite]{writeJuicer}} for details.
#' 
#' @param file Filename can be a local path or remote path. The remote path full list can be obtained from \url{http://aidenlab.org/data.html}.
#' @param verbose TRUE or FALSE. Whether print information.
#' 
#' @return A list object includes hic file information.
#' \item{\code{genomId}}{The genome id of the current hic file.}
#' \item{\code{resolution}}{The list of current hic file available resolution.}
#' \item{\code{pairs}}{The list of current hic file available pair.}
#' \item{\code{chromosomeSizes}}{A dataframe of chromosome size (chromosome, size). Can be used in juicer pre function for different genome.}
#' 
#' @examples 
#' library(FreeHiCLite)
#' 
#' ## Remote file location. The reomte file include downloading, it may take a while
#' remoteFilePath = 'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic'
#' 
#' ## Local file location
#' localFilePath = system.file('extdata', 'example.hic', package = 'FreeHiCLite')
#' 
#' juicerInfo <- readJuicerInformation(localFilePath)
#' print(str(juicerInfo))
#' 
#' \donttest{
#' juicerInfo <- readJuicerInformation(remoteFilePath)
#' print(str(juicerInfo))
#' }
#' 
#' @export
readJuicerInformation <- function(file, verbose = FALSE) {
    isHttp <- startsWith(file, "http")
    if (!isHttp) {
        file = normalizePath(file)
        stopifnot(file.exists(file))
    }
    # if (verbose) message('You are trying to access file from: ', file, '. Make sure you provide a
    # valid path.')
    
    # if (isHttp & verbose) message('Remote file can take a while.\nThe full list of remote files can
    # be found in http://aidenlab.org/data.html')
    ans <- hicDataInformation(file, isHttp)
    
    if (verbose) {
        message("File: ", file)
        message("GenomId: ", ans[["genomId"]])
        message("Available pair: ", paste0(ans[["pairs"]], collapse = ", "), ".")
        message("Hi-C resolution: ")
        message("    BP: ", paste0(ans[["resolution"]][["BP"]], collapse = ", "), ".")
        if (length(ans[["resolution"]][["FRAG"]]) > 0) {
            message("    FRAG: ", paste0(ans[["resolution"]][["FRAG"]], collapse = ", "), ".")
        } else {
            message("    FRAG not available.")
        }
    }
    
    class(ans) <- c("juicer", "information")
    return(ans)
}


#' @title 
#' 
#' Read fragment sites from the given hic file.
#' 
#' @description 
#' Read \href{https://github.com/aidenlab/juicer/wiki/Data}{.hic} file generated by \href{https://github.com/aidenlab/juicebox}{Juicebox}. Currently can be used for both local and remote file.
#' 
#' @details 
#' This function is heavily adopted from both java version juicebox \href{https://github.com/aidenlab/juicer/wiki/Data-Extraction}{Dump} function 
#' and c++ version \href{https://github.com/aidenlab/straw}{straw}.
#' 
#' @param file Filename can be a local path or remote path. The remote path full list can be obtained from \url{http://aidenlab.org/data.html}.
#' @param chromosomes A vector contains all the chromosomes. 
#' For example c('chr1', 'chr2')
#' @param verbose A logical value, indicate whether print message
#' 
#' @return A list contains all the chromosomes and their fragment sites.
#' 
#' @details 
#' 
#' The \href{https://github.com/aidenlab/juicer/wiki/Pre#restriction-site-file-format}{restriction site file} 
#' is required for juicer to provide a hic file with fragment level resolution.
#' 
#' 
#' @examples 
#' library(FreeHiCLite)
#' 
#' ## Remote file location. The reomte file include downloading, it may take a while
#' remoteFilePath = 'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic'
#' 
#' ## Local file location
#' localFilePath = system.file('extdata', 'example.hic', package = 'FreeHiCLite')
#' \donttest{
#' fragSites <- readJuicerFragmentSites(remoteFilePath, c('chr1', 'chr2'))
#' str(fragSites)
#' }
#' 
#' @export
readJuicerFragmentSites <- function(file, chromosomes, verbose = TRUE) {
    isHttp <- startsWith(file, "http")
    if (!isHttp) {
        file = normalizePath(file)
        stopifnot(file.exists(file))
    }
    if (verbose) 
        message("You are trying to access file from: ", file, ". Make sure you provide a valid path.")
    
    if (isHttp & verbose) 
        message("Remote file can take a while.\nThe full list of remote files can be found in http://aidenlab.org/data.html")
    
    ## process chromosome
    
    chromosomesChar <- NULL
    if (!is.null(chromosomes)) {
        chromosomesChar <- .chromosomes_clean(chromosomes)
    }
    return(hicDataFragSites(fileName = file, isHttp = isHttp, chromosomes = chromosomesChar))
}

.chromosomes_clean <- function(chromosomes) {
    chromosomesChar <- as.character(chromosomes)
    ans <- rep(NA, length(chromosomesChar))
    for (i in seq_along(chromosomesChar)) {
        ans[i] <- gsub("chr", "", tolower(trimws(chromosomesChar[i])))
    }
    sort(ans)
}


.check_resolution <- function(unit, resolution) {
    default <- list(BP = c(2500000, 1e+06, 5e+05, 250000, 1e+05, 50000, 25000, 10000, 5000), FRAG = c(500, 
        250, 100, 50, 20, 5, 2, 1))
    find <- resolution %in% default[[unit]]
    return(find)
}

.check_chr <- function(pairs, chromosomes) {
    for (pair in pairs) {
        tmp = unlist(strsplit(pair, "_"))
        if (!tmp[1] %in% chromosomes) {
            stop(tmp[1], " is not in chromosomes: ", paste(chromosomes, ", "))
        }
        if (!tmp[2] %in% chromosomes) {
            stop(tmp[2], " is not in chromosomes: ", paste(chromosomes, ", "))
        }
    }
    return(TRUE)
}

.parse_pair <- function(pairs) {
    ans <- rep(NA, length(pairs))
    for (i in seq_along(pairs)) {
        tmp <- unlist(strsplit(pairs[i], "_"))
        if (length(tmp) != 2) 
            return(FALSE)
        tmp <- paste0(.chromosomes_clean(tmp), collapse = "_")
        ans[i] = tmp
    }
    sort(unique(ans))
}

.chrosomes_to_pair <- function(chromosomes) {
    chromosomes = sort(unique(chromosomes))
    pair <- c()
    N <- length(chromosomes)
    for (i in 1:N) {
        for (j in i:N) {
            pair <- c(pair, paste0(chromosomes[i], "_", chromosomes[j]))
        }
    }
    pair
}

