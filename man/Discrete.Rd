\name{Discrete}
\alias{ddiscrete}
\alias{pdiscrete}
\alias{qdiscrete}
\alias{rdiscrete}
\title{Discrete Distribution}
\description{
  These functions provide information about the discrete distribution
  where the probability of the elements of \code{values} is proportional
  to the values given in \code{probs}, which are normalized to sum up to
  1.  \code{ddiscrete} gives the density, \code{pdiscrete} gives the
  distribution function, \code{qdiscrete} gives the quantile function
  and \code{rdiscrete} generates random deviates.
}
\usage{
ddiscrete(x, probs, values = 1:length(probs))
pdiscrete(q, probs, values = 1:length(probs))
qdiscrete(p, probs, values = 1:length(probs))
rdiscrete(n, probs, values = 1:length(probs), ...)
}
\arguments{
  \item{x,q}{vector or array of quantiles.}
  \item{p}{vector or array of probabilities.}
  \item{n}{number of observations.}
  \item{probs}{probabilities of the distribution.}
  \item{values}{values of the distribution.}
  \item{...}{ignored (only there for backwards compatibility)}
}
\details{
  The random number generator is simply a wrapper for
  \code{\link{sample}} and provided for backwards compatibility only.
}
\author{Andreas Weingessel and Friedrich Leisch}
\examples{
## a vector of length 30 whose elements are 1 with probability 0.2
## and 2 with probability 0.8.
rdiscrete (30, c(0.2, 0.8))

## a vector of length 100 whose elements are A, B, C, D.
## The probabilities of the four values have the relation 1:2:3:3
rdiscrete (100, c(1,2,3,3), c("A","B","C","D"))

}
\keyword{distribution}
