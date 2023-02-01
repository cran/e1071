\name{countpattern}
\alias{countpattern}
\title{Count Binary Patterns}
\usage{
countpattern(x, matching=FALSE)
}
\arguments{
 \item{x}{A matrix of binary observations}
 \item{matching}{If TRUE an additional vector is returned which stores
     which row belongs to which pattern}
}
\description{
Every row of the binary matrix \code{x} is transformed into a binary
pattern and these patterns are counted.
}
\value{
A vector of length \code{2\^ncol(x)} giving the number of times each
pattern occurs in the rows of \code{x}. The names of this vector are the
binary patterns. They are sorted according to their numeric value.
If \code{matching} is \code{TRUE}, a list of the following two vectors is returned.
  \item{pat}{Numbers of patterns as described above.}
  \item{matching}{Vector giving the position of the pattern of each row
      of \code{x} in \code{pat}.}
}

\author{Andreas Weingessel}

\examples{
xx <- rbind(c(1,0,0),c(1,0,0),c(1,0,1),c(0,1,1),c(0,1,1))
countpattern(xx)
countpattern(xx, matching=TRUE)
}

\keyword{multivariate}
