\name{kurtosis}
\alias{kurtosis}
\title{Kurtosis}
\description{
  Computes the kurtosis.
}
\usage{
kurtosis(x, na.rm = FALSE, type = 3)
}
\arguments{
  \item{x}{a numeric vector containing the values whose kurtosis is to
    be computed.}
  \item{na.rm}{a logical value indicating whether \code{NA} values
    should be stripped before the computation proceeds.}
  \item{type}{an integer between 1 and 3 selecting one of the algorithms
    for computing skewness detailed below.}
}
\details{
  If \code{x} contains missings and these are not removed, the skewness
  is \code{NA}.

  Otherwise, write \eqn{x_i} for the non-missing elements of \code{x},
  \eqn{n} for their number, \eqn{\mu}{mu} for their mean, \eqn{s} for
  their standard deviation, and
  \eqn{m_r = \sum_i (x_i - \mu)^r / n}{m_r = \sum_i (x_i - mu)^r / n}
  for the sample moments of order \eqn{r}.

  Joanes and Gill (1998) discuss three methods for estimating kurtosis:

  \describe{
    \item{Type 1:}{
      \eqn{g_2 = m_4 / m_2^2 - 3}.
      This is the typical definition used in many older textbooks.}
    \item{Type 2:}{
      \eqn{G_2 = ((n+1) g_2 + 6) * (n-1) / ((n-2)(n-3))}.
      Used in SAS and SPSS.
    }
    \item{Type 3:}{
      \eqn{b_2 = m_4 / s^4 - 3 = (g_2 + 3) (1 - 1/n)^2 - 3}.
      Used in MINITAB and BMDP.}
  }
  
  Only \eqn{G_2} (corresponding to \code{type = 2}) is unbiased under
  normality.
}
\value{
  The estimated kurtosis of \code{x}.
}
\references{
  D. N. Joanes and C. A. Gill (1998),
  Comparing measures of sample skewness and kurtosis.
  \emph{The Statistician}, \bold{47}, 183--189.
}
\examples{
x <- rnorm(100)
kurtosis(x)
}
\keyword{univar}
