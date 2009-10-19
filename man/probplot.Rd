\name{probplot}
\alias{probplot}
\alias{lines.probplot}
\title{Probability Plot}
\description{
  Generates a probability plot for a specified theoretical
  distribution, i.e., basically a
  \code{\link[stats]{qqplot}} where the y-axis is labeled with
  probabilities instead
  of quantiles. The function
  is mainly intended for teaching the concept of quantile plots.
}
\usage{
probplot(x, qdist=qnorm, probs=NULL, line=TRUE,
         xlab=NULL, ylab="Probability in \%", ...)
\S3method{lines}{probplot}(x, h=NULL, v=NULL, bend=FALSE, ...)
}
\arguments{
  \item{x}{A data vector for \code{probplot}, an object of class
    \code{probplot} for the \code{lines} method.}
  \item{qdist}{A character string or a function for the quantiles of the
    target distribution.}
  \item{probs}{Vector of probabilities at which horizontal lines should
    be drawn.}
  \item{line}{Add a line passing through the quartiles to the plot?}
  \item{xlab, ylab}{Graphical parameters.}
  \item{h}{The y-value for a horizontal line.}
  \item{v}{The x-value for a vertical line.}
  \item{bend}{If \code{TRUE}, lines are ``bent'' at the quartile line, else
    regular \code{abline}s are added. See examples.}
  \item{\dots}{Further arguments for \code{qdist} and graphical
    parameters for lines.}
}

\author{Friedrich Leisch}
\seealso{\code{\link[stats]{qqplot}}}
\examples{
## a simple example
x <- rnorm(100, mean=5)
probplot(x)

## the same with horizontal tickmarks at the y-axis
opar <- par("las")
par(las=1)
probplot(x)

## this should show the lack of fit at the tails
probplot(x, "qunif")

## for increasing degrees of freedom the t-distribution converges to
## normal
probplot(x, qt, df=1)
probplot(x, qt, df=3)
probplot(x, qt, df=10)
probplot(x, qt, df=100)

## manually add the line through the quartiles
p <- probplot(x, line=FALSE)
lines(p, col="green", lty=2, lwd=2)


## Make the line at prob=0.5 red
lines(p, h=0.5, col="red")

### The following use the estimted distribution given by the green
### line:

## What is the probability that x is smaller than 7?
lines(p, v=7, bend=TRUE, col="blue")

## Median and 90\% confidence interval
lines(p, h=.5, col="red", lwd=3, bend=TRUE)
lines(p, h=c(.05, .95), col="red", lwd=2, lty=3, bend=TRUE)

par(opar)
}
\keyword{hplot}
