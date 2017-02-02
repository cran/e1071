\name{rbridge}
\alias{rbridge}
\title{Simulation of Brownian Bridge}
\usage{
rbridge(end = 1, frequency = 1000)
}
\arguments{
\item{end}{the time of the last observation.}
\item{frequency}{the number of observations per unit of time.}
}
\description{
\code{rwiener} returns a time series containing a simulated realization
of the Brownian bridge on the interval [0,\code{end}]. If W(t) is a
Wiener process, then the Brownian bridge is defined as W(t) - t W(1).
}
\seealso{
rwiener
}
\examples{
# simulate a Brownian bridge on [0,1] and plot it

x <- rbridge()
plot(x,type="l")
}
\keyword{distribution}
