\name{hamming.window}

\title{Computes the Coefficients of a Hamming Window.}

\usage{hamming.window(n)}

\alias{hamming.window}


\arguments{
\item{n}{The length of the window.}
}

\description{The filter coefficients \eqn{w_i}{w(i)} of a Hamming
window of length 
\code{n} are computed according to the formula
\deqn{w_i = 0.54 - 0.46 \cos\frac{2\pi i}{n-1}}{
w(i) = 0.54 - 0.46*cos(2*pi*i/(n-1))}
}

\value{A vector containing the filter coefficients.}

\references{For a definition of the Hamming window, see for example\cr
Alan V. Oppenheim and Roland W. Schafer: "Discrete-Time Signal
Processing", Prentice-Hall, 1989.}

\author{Andreas Weingessel}

\seealso{stft, hanning.window}

\examples{hamming.window(10)

x<-rnorm(500)
y<-stft(x, wtype="hamming.window")
plot(y)
}

\keyword{ts}
