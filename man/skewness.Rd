\name{skewness}
\alias{skewness}
\title{Skewness}
\usage{
skewness(x, na.rm=FALSE)
}
\arguments{
\item{x}{a numeric vector containing the values whose skewness is to
be computed.}
\item{na.rm}{a logical value indicating whether NA values should be
stripped before the computation proceeds.}
}
\description{
Computes the skewness of the values of \code{x}.

If \eqn{N = length(x)}, then the skewness of \eqn{x} is defined as

\deqn{N^(-1) sd(x)^(-3) sum_i (x_i - mean(x))^3}  
}
\value{
Returns the skewness of x.
}
\examples{
x<-rnorm(100)
skewness(x)
}
\keyword{univar}

