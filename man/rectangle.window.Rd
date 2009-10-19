\name{rectangle.window}
\title{Computes the Coefficients of a Rectangle Window.}

\usage{rectangle.window(n)}

\alias{rectangle.window}


\arguments{
\item{n}{The length of the window.}
}

\description{Returns the filter coefficients of a rectangle
window. That is a vector of \code{n} 1.

The purpose of this function is just to have a name for the R command
\code{rep (1, n)}. 
}

\value{A vector of length \code{n} filled with 1.}

\author{Andreas Weingessel}

\seealso{stft}

\examples{x<-rnorm(500)
y<-stft(x, wtype="rectangle.window")
plot(y)
}

\keyword{ts}
