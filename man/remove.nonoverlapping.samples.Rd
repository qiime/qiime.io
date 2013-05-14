\name{remove.nonoverlapping.samples}
\alias{remove.nonoverlapping.samples}
\title{Ensure tables imported from QIIME contain same samples in the same order}
\usage{
  remove.nonoverlapping.samples(map = NULL, otus = NULL,
    taxa = NULL, distmat = NULL)
}
\arguments{
  \item{map}{object}

  \item{otus}{object}

  \item{taxa}{object}

  \item{distmat}{object}
}
\value{
  A list with your \code{\link{data.frame}} or
  \code{\link{matrix}} objects, with samples ordered and
  intersected.
}
\description{
  Some details about which files. Links to
  \code{\link{load.qiime.distance.matrix}}
  \code{\link{load.qiime.taxon.table}}
  \code{\link{load.qiime.mapping.file}}
  \code{\link{load.qiime.otu.table}}
}
\examples{
# Some R code examples here.
# Need example files to include in package.
}
\references{
  \url{http://qiime.org/}

  ``QIIME allows analysis of high-throughput community
  sequencing data.'' Nature Methods, 2010;
  doi:10.1038/nmeth.f.303
}
\seealso{
  \code{\link{load.qiime.distance.matrix}}

  \code{\link{load.qiime.taxon.table}}

  \code{\link{load.qiime.mapping.file}}

  \code{\link{load.qiime.otu.table}}
}

