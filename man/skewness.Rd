\name{skewness}
\alias{skewness}
\title{Skewness}
\description{
  Computes the skewness.
}
\usage{
skewness(x, na.rm=FALSE)
}
\arguments{
  \item{x}{a numeric vector containing the values whose skewness is to
    be computed.}
  \item{na.rm}{a logical value indicating whether \code{NA} values
    should be stripped before the computation proceeds.}
}
\details{
  If \eqn{N = \mathrm{length}(x)}{N = length(x)}, then the skewness of
  \eqn{x} is defined as 
  \deqn{N^{-1} \mathrm{sd}(x)^{-3} \sum_i (x_i - \mathrm{mean}(x))^3.}{
    N^(-1) sd(x)^(-3) sum_i (x_i - mean(x))^3.}
}
\value{
  The skewness of \code{x}.
}
\examples{
x <- rnorm(100)
skewness(x)
}
\keyword{univar}
