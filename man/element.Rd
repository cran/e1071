\name{element}
\title{Extract Elements of an Array}
\usage{
element(x, i)
}
\alias{element}
\arguments{
 \item{x}{Array of arbitrary dimensionality.}
 \item{i}{Vector of the same length as \code{x} has dimension.}
}
\description{
  Returns the element of \code{x} specified by \code{i}.
}
\author{Friedrich Leisch}
\seealso{\link{Extract}}
\examples{
x <- array(1:20, dim=c(2,5,2))
element(x, c(1,4,2))
}
\keyword{array}
