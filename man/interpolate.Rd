\name{interpolate}
\title{Interpolate Values of Array}
\usage{
interpolate(x, a, adims=lapply(dimnames(a), as.numeric),
            method="linear")
}
\alias{interpolate}
\arguments{
 \item{x}{Matrix of values at which interpolation shall take place.}
 \item{a}{Array of arbitrary dimension.}
 \item{adims}{List of the same structure as \code{dimnames(a)}.}
 \item{method}{Interpolation method, one of \code{"linear"} or
   \code{"constant"}.}
}
\description{
  For each row in matrix \code{x}, the hypercube of \code{a} containing
  this point is searched. The corners of the hypercube are
  linearly interpolated. By default, \code{dimnames(a)} is taken to
  contain the coordinate values for each point in \code{a}. This can be
  overridden using \code{adims}. If \code{method=="constant"}, the value
  of the ``lower left'' corner of the hypercube is returned.
}
\author{Friedrich Leisch}
\seealso{\code{\link{approx}}, \code{\link{spline}}}
\examples{
x <- seq(0,3,0.2)
z <- outer(x,x, function(x,y) sin(x*y))
dimnames(z) <- list(x,x)
sin(1.1*2.1)
interpolate(c(1.1, 2.1),z)
}
\keyword{arith}
\keyword{multivariate}

