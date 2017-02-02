\name{moment}
\alias{moment}
\title{Statistical Moment}
\description{
  Computes the (optionally centered and/or absolute) sample moment of a
  certain order.
}
\usage{
moment(x, order=1, center=FALSE, absolute=FALSE, na.rm=FALSE)
}
\arguments{
 \item{x}{a numeric vector containing the values whose moment is to be
   computed.}
 \item{order}{order of the moment to be computed, the default is to
   compute the first moment, i.e., the mean.}
 \item{center}{a logical value indicating whether centered moments are
   to be computed.}
 \item{absolute}{a logical value indicating whether absolute moments are
   to be computed.}
 \item{na.rm}{a logical value indicating whether \code{NA} values
   should be stripped before the computation proceeds.}
}
\details{
  When \code{center} and \code{absolute} are both \code{FALSE}, the
  moment is simply \code{sum(x ^ order) / length(x)}.
}
\author{Kurt Hornik and Friedrich Leisch}
\seealso{
  \code{\link{mean}},
  \code{\link{var}}}
\examples{
x <- rnorm(100)

## Compute the mean
moment(x)
## Compute the 2nd centered moment (!= var)
moment(x, order=2, center=TRUE)

## Compute the 3rd absolute centered moment
moment(x, order=3, center=TRUE, absolute=TRUE)
}
\keyword{univar}
