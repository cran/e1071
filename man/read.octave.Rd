\name{read.octave}
\alias{read.octave}
\title{Read Octave Data File}
\usage{
read.octave (file, quiet=FALSE)
}
\arguments{
\item{file}{the name of a Octave data file to read.}
\item{quiet}{do not print information on type and size of the object
read in.}
}
\description{
Read a vector or matrix from an Octave ASCII data file created using
the command \code{save -ascii} in Octave. \code{read.octave} reads only
the first object found in the file and ignores the rest.
}
\value{
\code{read.octave} returns an object of type vector or matrix.
}
\seealso{
scan
read.table
}
\keyword{interface}

