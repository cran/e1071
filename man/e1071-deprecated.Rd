\name{e1071-deprecated}
\alias{e1071-deprecated}
\alias{read.octave}
\title{Deprecated Functions in Package e1071}
\description{
  These functions are provided for compatibility with older versions of
  package \pkg{e1071} only, and may be defunct as soon as of the next
  release.
}
\usage{
read.octave(file, quiet = FALSE)
}
\arguments{
  \item{file}{a character string giving the name of a file in Octave
    text data format (as created by \code{save -ascii}) to read.} 
  \item{quiet}{do not print information on type and size of the object
    read in.}
}
\details{
  \code{read.octave} reads only the first object found in the file and
  ignores the rest.
}
\value{
  \code{read.octave} returns an object of type vector or matrix.
}
\seealso{
  \code{\link{Deprecated}};
  \code{\link[foreign]{read.octave}}
}
\keyword{misc}
