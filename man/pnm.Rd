\name{pnm}
\title{Portable Anymap Images}
\usage{
read.pnm(file)
plot(pnmobj, xlab="", ylab="", axes=FALSE, ...)
write.pgm(pnmobj, file="Rimage.pgm", forceplain=FALSE)
}
\alias{write.pgm}
\alias{read.pnm}
\alias{plot.pnm}
\arguments{
  \item{pnmobj}{an object of class \code{pnm}}
  \item{file}{name of the pnm file}
  \item{forceplain}{If \code{TRUE}, an ASCII pgm file is written. Default
    is to write a binary (raw) pgm file.}
}
\description{
  \code{read.pnm} reads a pnm file and loads the image into an
  object of class \code{pnm}, which is basically a
  matrix (pbm, pgm) or a 3-dimensional arry with 3 matrices for red,
  green and blue.

  \code{write.pgm} writes an object of class \code{pnm} to a
  pgm file. Generalizations for writing pbm and ppm files are yet to be
  written.

  \code{plot.pnm} plots a pnm object using the command
  \code{\link{image}}. The only difference is that the element \code{[1,1]} of
  \code{pnmobj} is plotted as the upper left corner (plain
  \code{\link{image}} would plot \code{[1,1]} as the lower left corner.
}
\value{
  \code{read.pnm} returns an object of class \code{pnm}. The maximum
  value (of either grey or the red/green/blue channels) is stored as
  attribute \code{"maxval"}). The type (pbm, pgm or pnm) is stred as
  attribute \code{"type"}).
}
\author{Friedrich Leisch}
\seealso{\code{\link{image}}}
\keyword{file,hplot}
