\name{kurtosis}
\alias{kurtosis}
\title{Kurtosis}
\usage{
kurtosis(x, na.rm=FALSE)
}
\arguments{
\item{x}{a numeric vector containing the values whose kurtosis is to
be computed.}
\item{na.rm}{a logical value indicating whether NA values should be
stripped before the computation proceeds.}
}
\description{
Computes the kurtosis of the values of \code{x}.

If \eqn{N = length(x)}, then the kurtosis of \eqn{x} is defined as
\deqn{N^(-1) var(x)^(-2) sum_i (x_i - mean(x))^4 - 3}  
}
\value{
Returns the kurtosis of x.
}
\examples{
x<-rnorm(100)
kurtosis(x)
}
\keyword{univar}
