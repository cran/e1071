\name{kurtosis}
\alias{kurtosis}
\title{Kurtosis}
\description{
  Computes the kurtosis.
}
\usage{
kurtosis(x, na.rm=FALSE)
}
\arguments{
  \item{x}{a numeric vector containing the values whose kurtosis is to
    be computed.}
  \item{na.rm}{a logical value indicating whether \code{NA} values
    should be stripped before the computation proceeds.}
}
\details{
  If \eqn{N = \mathrm{length}(x)}{N = length(x)}, then the kurtosis of
  \eqn{x} is defined as
  \deqn{N^{-1} \mathrm{sd}(x)^{-4}
    \sum_i(x_i - \mathrm{mean}(x))^4 - 3.}{
    N^(-1) sd(x)^(-4) sum_i (x_i - mean(x))^4 - 3.}
}
\value{
  The kurtosis of \code{x}.
}
\examples{
x <- rnorm(100)
kurtosis(x)
}
\keyword{univar}
