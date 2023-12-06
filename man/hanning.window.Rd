\name{hanning.window}

\title{Computes the Coefficients of a Hanning Window.}

\usage{hanning.window(n)}

\alias{hanning.window}


\arguments{
\item{n}{The length of the window.}
}

\description{The filter coefficients \eqn{w_i}{w(i)} of a Hanning
window of length 
\code{n} are computed according to the formula
\deqn{w_i = 0.5 - 0.5 \cos\frac{2\pi i}{n-1}}{
w(i) = 0.5 - 0.5*cos(2*pi*i/(n-1))}
}

\value{A vector containing the filter coefficients.}

\references{For a definition of the Hanning window, see for example\cr
Alan V. Oppenheim and Roland W. Schafer: "Discrete-Time Signal
Processing", Prentice-Hall, 1989.}

\author{Andreas Weingessel}

\seealso{stft, hamming.window}

\examples{hanning.window(10)

x<-rnorm(500)
y<-stft(x, wtype="hanning.window")
plot(y)
}

\keyword{ts}
