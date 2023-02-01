\name{stft}
\title{Computes the Short Time Fourier Transform of a Vector}

\usage{stft(X, win=min(80,floor(length(X)/10)), inc=min(24,
floor(length(X)/30)), coef=64, wtype="hanning.window")}  
\alias{stft}

\arguments{
\item{X}{The vector from which the stft is computed.}
\item{win}{Length of the window. For long vectors the default window
size is 80, for short vectors the window size is chosen so that 10
windows fit in the vector.}
\item{inc}{Increment by which the window is shifted. For long vectors
the default increment is 24, for short vectors the increment is chosen
so that 30 increments fit in the vector.} 
\item{coef}{Number of Fourier coefficients}
\item{wtype}{Type of window used}
}

\description{This function computes the Short Time Fourier Transform
of a given vector \code{X}. 

First, time-slices of length \code{win} are
extracted from the vector. The shift of one time-slice to the next one is
given by \code{inc}. The values of these time-slices are smoothed by
mulitplying them with a window function specified in \code{wtype}. For
the thus obtained windows, the Fast Fourier Transform is computed.}

\value{Object of type stft. Contains the values of the stft and
information about the parameters. 

\item{values}{A matrix containing the results of the stft. Each row of
the matrix contains the \code{coef} Fourier coefficients of one
window.} 
\item{windowsize}{The value of the parameter \code{win}}
\item{increment}{The value of the parameter \code{inc}}
\item{windowtype}{The value of the parameter \code{wtype}}
}

\author{Andreas Weingessel}

\seealso{plot.stft}

\examples{x<-rnorm(500)
y<-stft(x)
plot(y)
}

\keyword{ts}

