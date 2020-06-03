#' @title
#' 
#' Write to short format with score
#' 
#' @description 
#' 
#' This function will write contact matrixes into \href{https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format}{short with score format}.
#' 
#' 
#' @param contactMatrixsList A list contains all the contact matrix or a single matrix. The list name should be pair as form '1_1' or '1_2'
#' @param pairs Which pairs of the chromosomes. If \code{contactMatrixList} is a list, with name, this argument can be NULL. 
#' If \code{contactMatrixList} is a matrix, it must has pair which format "A_B".
#' @param file File name to store the contact matrix
#' @param overwrite If file exists, whether overwrite the file or append file. If TRUE, the function will remove the file and 
#' re-write it. If FALSE, function will append the result to the existing file.
#' 
#' @details 
#' 
#' The \href{https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format}{short with score format} contains 9 columns.
#' 
#' \enumerate{
#' \item str1 = strand (0 for forward, anything else for reverse)
#' \item chr1 = the first chromosome (must be a chromosome in the genome)
#' \item pos1 = the first position 
#' \item frag1 = restriction site fragment
#' \item str2 = strand (0 for forward, anything else for reverse)
#' \item chr2 = the second chromosome (must be a chromosome in the genome)
#' \item pos2 = the second position 
#' \item frag2 = restriction site fragment
#' \item score = the score imputed to this read
#' }
#' 
#' If not using the restriction site file option \code{\link[FreeHiCLite]{readJuicerFragmentSites}}, both frag1 and frag2 will be ignored. By default, the function gives them 0.
#' 
#' \code{java -jar juicer_tool.jar -r 5000 -n <infile> <outfile> <genomeID>}
#' 
#' If \code{<genomeID>} not in the list of hg18, hg19, hg38, dMel, mm9, mm10, 
#' anasPlat1, bTaurus3, canFam3, equCab2, galGal4, Pf3D7, sacCer3, sCerS288c, susScr3, or TAIR10. 
#' You must provide a chrom.sizes file path that lists on each line the name and size of the chromosomes.
#' 
#' Use \code{\link[FreeHiCLite]{readJuicerInformation}} with the given hic file, you can extract the chromosome size. 
#' The resulting list contains a list \code{chromosomeSizes}.
#' 
#' 
#' @return write a file into disk
#' 
#' @examples 
#' 
#' # From existing object
#' 
#' library(FreeHiCLite)
#'  
#' ## Local file location
#' localFilePath = system.file('extdata', 'example.hic', package = 'FreeHiCLite')
#' 
#' ## Pairs needs to be extracted
#' pairs = c('1_1', '1_2')
#' unit = 'BP'
#' resolution = 500000L
#' 
#' ## pass chrosomes into function, it will contains all the interaction pairs
#' 
#' dat <- readJuicer(file=localFilePath, chromosomes=NULL, 
#' pairs = pairs, unit=unit, resolution=resolution)
#' 
#' contactMatrixsList <- dat[['contact']]
#' 
#' \dontshow{.old_wd <- setwd(tempdir())}
#' writeJuicer(contactMatrixsList, file = 'test.txt')
#' \dontshow{setwd(.old_wd)}
#' 
#' 
#' @importFrom utils menu
#' @importFrom methods is
#' @export
writeJuicer <- function(contactMatrixsList, pairs = NULL, file, overwrite = TRUE) {
    if (is.null(pairs)) pairs <- names(contactMatrixsList)
    
    if (file.exists(file)) {
        if (overwrite) 
            file.remove(file) else {
            tmp = utils::menu(c("Overwrite", "Append"), title = paste(file, "exists."))
            if (tmp == 1) 
                file.remove(file)
        }
    }
    
    if (methods::is(contactMatrixsList, "matrix")) {
        stopifnot(!is.null(pairs))
        .write.juicer(contactMatrixsList, pairs, file, append = TRUE)
    } else {
        for (pair in pairs) {
            .write.juicer(contactMatrixsList[[pair]], pair, file, append = TRUE)
        }
    }
}


#' @title
#' 
#' Convert contact matrix
#' 
#' @description 
#' 
#' This function can covert a contact matrix (x, y, counts) format into a \href{https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format}{short with score format}.
#' 
#' @param contactMatrix A contact matrix, the first column is the first chromosome location, and second column is the second chromosome location, the third column is count number
#' @param chr1 The first chromosome name. For example, 1 or 'chr1'.
#' @param chr2 The second chromosome name. For example, 1 or 'chr1'
#' 
#' @return A data.frame has 9 columns, and it can be used directly by juicer \href{https://github.com/aidenlab/juicer/wiki/Pre}{Pre}.
#' 
#' 
#' @seealso \code{\link[FreeHiCLite]{writeJuicer}}.
#' 
#' @examples
#' 
#' ## Create a random contact matrix
#' 
#' contactMatrix <- structure(
#' c(0L, 0L, 50L, 0L, 50L, 100L, 0L, 50L, 50L, 100L, 100L, 
#' 100L, 526L, 123L, 499L, 36L, 213L, 562L), .Dim = c(6L, 3L), 
#' .Dimnames = list(NULL, c('x', 'y', 'counts')))
#' 
#' contactMatrix
#' 
#' ## chromosomes names
#' 
#' chr1 = 'chr1'
#' chr2 = 'chr2'
#' 
#' ## Convert to hic
#' (df = convertJuicer(contactMatrix, chr1, chr2))
#' 
#' @export
convertJuicer <- function(contactMatrix, chr1, chr2) {
    df <- data.frame(str1 = 0, chr1 = chr1, pos1 = contactMatrix[, 1], frag1 = 0, str2 = 0, chr2 = chr2, 
        pos2 = contactMatrix[, 2], frag2 = 1, score = contactMatrix[, 3])
    return(df)
}


#' @importFrom utils write.table
.write.juicer <- function(contactMatrix, pair, file, append = FALSE) {
    
    tmp <- unlist(strsplit(pair, "_"))
    chr1 <- .chromosomes_clean(tmp[1])
    chr2 <- .chromosomes_clean(tmp[2])
    
    mat <- convertJuicer(contactMatrix, chr1, chr2)
    
    utils::write.table(mat, file = file, append = append, quote = FALSE, sep = "\t", row.names = FALSE, 
        col.names = FALSE)
}
