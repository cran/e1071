\name{rwiener}
\alias{rwiener}
\title{Simulation of Wiener Process}
\usage{
rwiener(end = 1, frequency = 1000)
}
\arguments{
\item{end}{the time of the last observation.}
\item{frequency}{the number of observations per unit of time.}
}
\description{
\code{rwiener} returns a time series containing a simulated realization
of the Wiener process on the interval [0,\code{end}]
}
\examples{
# simulate a Wiener process on [0,1] and plot it

x <- rwiener()
plot(x,type="l")
}
\keyword{distribution}
