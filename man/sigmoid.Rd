\name{sigmoid}
\title{The logistic function and derivatives}
\alias{sigmoid}
\alias{dsigmoid}
\alias{d2sigmoid}
\usage{
sigmoid(x)
dsigmoid(x)
d2sigmoid(x)
}
\description{Sigmoid \eqn{1/(1 + exp(-x))}, first and second derivative.}
\author{Friedrich Leisch}
\examples{
plot(sigmoid, -5, 5, ylim=c(-.2, 1))
plot(dsigmoid, -5, 5, add=T, col=2)
plot(d2sigmoid, -5, 5, add=T, col=3)
}
\keyword{math}
