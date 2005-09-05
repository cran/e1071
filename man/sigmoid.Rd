\name{sigmoid}
\alias{sigmoid}
\alias{dsigmoid}
\alias{d2sigmoid}
\title{The Logistic Function and Derivatives}
\description{
  Sigmoid \eqn{1/(1 + \exp(-x))}{1/(1 + exp(-x))}, first and second
  derivative.}
\usage{
sigmoid(x)
dsigmoid(x)
d2sigmoid(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\author{Friedrich Leisch}
\examples{
plot(sigmoid, -5, 5, ylim = c(-.2, 1))
plot(dsigmoid, -5, 5, add = TRUE, col = 2)
plot(d2sigmoid, -5, 5, add = TRUE, col = 3)
}
\keyword{math}

