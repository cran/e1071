\name{hamming.distance}
\alias{hamming.distance}
\title{Hamming Distances of Vectors}
\usage{ hamming.distance(x, y)
}
\arguments{
\item{x}{a vector or matrix.}
\item{y}{an optional vector.}
}
\description{ If both \code{x} and
\code{y} are vectors, \code{hamming.distance} returns the Hamming
distance (number of different bytes) between this two vectors. If
\code{x} is a matrix, the Hamming distances between the rows of \code{x}
are computed and \code{y} is ignored.
}
\examples{
x <- c(1, 0, 0)
y <- c(1, 0, 1)
hamming.distance(x, y)
z <- rbind(x,y)
rownames(z) <- c("Fred", "Tom")
hamming.distance(z)
}
\keyword{multivariate}
