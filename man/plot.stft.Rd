\name{plot.stft}
\title{Plots an Object of Type stft.}

\usage{plot.stft(Y, col=gray(63:0/63))}
\alias{plot.stft}

\arguments{
\item{Y}{An object of class \code{stft} as obtained by the function
\code{stft}.}  
\item{col}{An optional colormap. By default 64 gray values are used,
where white corresponds to the minimum value and black to the
maximum.}
}

\description{An object of class \code{stft} is plotted as a gray scale
image. The x-axis corresponds to time, the y-axis to frequency. If the
default colormap is used, dark regions in the plot correspond to high
values at the particular time/frequency location.}

\value{No return value. This function is only for plotting.}

\author{Andreas Weingessel}

\seealso{stft}

\examples{x<-rnorm(500)
y<-stft(x)
plot(y)
}

\keyword{ts}
