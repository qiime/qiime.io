\name{load.qiime.otu.table}
\alias{load.qiime.distance.matrix}
\alias{load.qiime.mapping.file}
\alias{load.qiime.otu.table}
\alias{load.qiime.taxon.table}
\title{Read files created by the QIIME pipeline.}
\usage{
  load.qiime.otu.table(filepath, include.lineages = FALSE)

  load.qiime.mapping.file(filepath)

  load.qiime.taxon.table(filepath)

  load.qiime.distance.matrix(filepath)
}
\arguments{
  \item{filepath}{(Required). The character string file
  path to the legacy QIIME OTU file that you want to
  import.  Must be a path suitable for
  \code{\link[base]{scan}} or
  \code{\link[utils]{read.table}}.  Can also be a
  \code{\link[base]{url}} or R \code{\link{connection}}.}

  \item{include.lineages}{(Optional). A logical with length
  one.  Should lineages be included in the returned object.
  Default is \code{FALSE}.  This is only relevant for
  \code{load.qiime.otu.table}}
}
\value{
  A \code{\link{data.frame}} or \code{\link{matrix}}
}
\description{
  Some details about qiime legacy files here.
}
\examples{
# mytoutaxfilepath <- "example/path/qiime/legacy/file"
 # w = read.qiime.table(mytoutaxfilepath)
 # x = read.qiime.table(mapfilepath)
 # y = read.qiime.table(taxfilepath)
 # z = read.qiime.table(distancefilepath)
 # dim(x)
 # head(x)[, 1:10]
 ## Other stuff
}
\references{
  \url{http://qiime.org/}

  ``QIIME allows analysis of high-throughput community
  sequencing data.'' Nature Methods, 2010;
  doi:10.1038/nmeth.f.303
}
\seealso{
  \code{\link{remove.nonoverlapping.samples}}
}

