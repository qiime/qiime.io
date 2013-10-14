###############################################################################
#' Read files created by the QIIME pipeline.
#'
#' Some details about qiime legacy files here.
#' 
#' @param filepath (Required). The character string file path to the legacy 
#'  QIIME OTU file that you want to import.
#'  Must be a path suitable for 
#'  \code{\link[base]{scan}} or \code{\link[utils]{read.table}}.
#'  Can also be a \code{\link[base]{url}}
#'  or R \code{\link{connection}}.
#'
#' @param include.lineages (Optional). A logical with length one.
#'  Should lineages be included in the returned object.
#'  Default is \code{FALSE}.
#'  This is only relevant for \code{load.qiime.otu.table}
#'
#' @return A \code{\link{data.frame}} or \code{\link{matrix}}
#'
#' @seealso
#'  \code{\link{remove.nonoverlapping.samples}}
#'
#' @references \url{http://qiime.org/}
#'
#' ``QIIME allows analysis of high-throughput community sequencing data.''
#' Nature Methods, 2010; doi:10.1038/nmeth.f.303
#'
#' @aliases load.qiime.distance.matrix
#' @aliases load.qiime.taxon.table
#' @aliases load.qiime.mapping.file
#' @aliases load.qiime.otu.table
#' 
#' @rdname input-methods
#' @export
#' @examples
#'  # mytoutaxfilepath <- "example/path/qiime/legacy/file"
#'  # w = read.qiime.table(mytoutaxfilepath)
#'  # x = read.qiime.table(mapfilepath)
#'  # y = read.qiime.table(taxfilepath)
#'  # z = read.qiime.table(distancefilepath)
#'  # dim(x)
#'  # head(x)[, 1:10]
#'  ## Other stuff
"load.qiime.otu.table" <- function(filepath, include.lineages=FALSE){
  
  otus <- read.qiime.table(filepath, as.data.frame=TRUE)
  # drop "Consensus Lineage" column if present
  if(otu.table.has.metadata(colnames(otus))){
    C <- ncol(otus)
    lineages <- as.character(otus[, C])
    otus <- otus[, -C]
  } else {
    lineages <- NULL
  }
  otus <- as.matrix(t(otus))
  
  if(include.lineages){
    return(list(otus=otus, lineages=lineages))
  } else {
    return(otus=otus)
  }
}

#' @rdname input-methods
#' @export
"load.qiime.mapping.file" <- function(filepath){
  return(read.qiime.table(filepath, as.data.frame=TRUE))
}

#' @rdname input-methods
#' @export
"load.qiime.taxon.table" <- function(filepath){
  taxa <- as.matrix(t(read.table(filepath, sep='\t', header=TRUE,
                                 row.names=1, check.names=FALSE, quote='"')))
  return(taxa)
}

#' @rdname input-methods
#' @export
"load.qiime.distance.matrix" <- function(filepath){
  d <- as.matrix(read.table(filepath,sep='\t', header=TRUE,
                            row.names=1, check.names=FALSE, quote='"'))
  return(d)
}

#
#' Ensure tables imported from QIIME contain same samples in the same order
#'
#' Some details about which files.
#' Links to
#' \code{\link{load.qiime.distance.matrix}}
#' \code{\link{load.qiime.taxon.table}}
#' \code{\link{load.qiime.mapping.file}}
#' \code{\link{load.qiime.otu.table}}
#'
#' @param map object
#' @param otus object
#' @param taxa object
#' @param distmat object
#' 
#' @return A list with your
#'  \code{\link{data.frame}} or \code{\link{matrix}} objects,
#'  with samples ordered and intersected.
#' 
#' @seealso
#'  \code{\link{load.qiime.distance.matrix}}
#'  
#'  \code{\link{load.qiime.taxon.table}}
#'  
#'  \code{\link{load.qiime.mapping.file}}
#'  
#'  \code{\link{load.qiime.otu.table}}
#'
#'
#' @references \url{http://qiime.org/}
#'
#' ``QIIME allows analysis of high-throughput community sequencing data.''
#' Nature Methods, 2010; doi:10.1038/nmeth.f.303
#'
#' @export
#' @examples
#' # Some R code examples here.
#' # Need example files to include in package.
"remove.nonoverlapping.samples" <- function(map=NULL, otus=NULL, taxa=NULL, distmat=NULL){
  IDs <- NULL
  objects <- list(map=map, otus=otus, taxa=taxa, distmat=distmat)
  
  # find overlapping samples in all tables
  for(obj in objects){
    if(!is.null(obj)) {
      if(is.null(IDs)){
        IDs <- rownames(obj)
      } else {
        IDs <- intersect(rownames(obj), IDs)
      }
    }
  }
  
  # drop non-overlapping samples 
  for(i in 1:length(objects)){
    if(!is.null(objects[[i]])) {
      objects[[i]] <- objects[[i]][IDs,,drop=F]
      # for mapping file, drop any empty levels from factors that might
      # have occurred due to dropped samples
      if(i == 1) objects[[i]] <- droplevels(objects[[i]])
      # for distance matrix, get subset of columns too
      if(i == 4) objects[[i]] <- objects[[i]][,IDs,drop=F]
    }
  }
  
  return(objects)
}
###############################################################################
# Internal functions not exported to namespace for users.
###############################################################################

# The core internal file-reading function
#' @keywords internal
"read.qiime.table" <- function(filepath, as.data.frame=FALSE){
  header.index <- get.header.index(filepath)
  # read the header
  f <- file(filepath, 'r')
  header <- scan(filepath, what='character', sep='\t', comment.char='',
                 skip=header.index-1, quote='"',
                 nlines=1, quiet=TRUE)
  # Close the connection
  close(f)
  # read the rest of the table
  datatable <- read.table(filepath, sep='\t',
                          skip=header.index, comment.char='#', quote='"',
                          header=FALSE, row.names=1, check.names=FALSE,
                          strip.white=TRUE)
  
  # set column names using header
  colnames(datatable) <- header[-1]
  
  if(!as.data.frame){
    datatable <- as.matrix(datatable)
  }
  return(datatable)
}

#
# TRUE if last column is "Consensus Lineage" or "OTU Metadata"
#' @keywords internal
"otu.table.has.metadata" <- function(headers){
  C <- length(headers)
  has.metadata <- grepl('consensus[ ]lineage|otu[ ]*metadata',
                        headers[C], ignore.case=TRUE)
  return(has.metadata)
}

# returns the index of the header line
# note: lines after the header may be comments with '#'
# read.table should therefore be called with (skip=header.index, comment.char='#')
#' @keywords internal
"get.header.index" <- function(filepath){
  ncolumns.per.line <- NULL
  # read lines until the first line without a '#'
  # for each line, obtain the number of tab-delimited columns
  linecount <- 0
  start.character <- '#'
  while(start.character == '#'){
    linecount <- linecount + 1
    f <- file(filepath,'r') # open file in read mode
    line <- scan(f,what='character',skip=linecount-1,nlines=1, sep='\t', quiet=TRUE)
    close(f)
    # ncolumns is the number of entries in this line
    # not including trailing empties
    ncolumns <- max(which(sapply(line,nchar) > 0))
    ncolumns.per.line <- c(ncolumns.per.line, ncolumns)
    start.character <- substring(line[1],1,1)
  }
  # first non-comment line gives the number of columns
  C <- ncolumns.per.line[linecount]
  if(linecount == 1){
    # if there are no comment lines, then the first line is the header
    header.index <- 1
  } else {
    if(any(ncolumns.per.line[-linecount] == C)){
      # if there is a comment line with the correct number of columns,
      # it is the header
      header.index <- max(which(ncolumns.per.line[-linecount] == C))
    } else {
      # if there is no comment line with the correct number of columns,
      # the first non-comment line is the header
      header.index <- linecount
    }
  }
  return(header.index)
}
###############################################################################
