\name{skewness}
\alias{skewness}
\title{Skewness}
\description{
  Computes the skewness.
}
\usage{
skewness(x, na.rm = FALSE, type = 3)
}
\arguments{
  \item{x}{a numeric vector containing the values whose skewness is to
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

  Joanes and Gill (1998) discuss three methods for estimating skewness:

  \describe{
    \item{Type 1:}{
      \eqn{g_1 = m_3 / m_2^{3/2}}{g_1 = m_3 / m_2^(3/2)}.
      This is the typical definition used in many older textbooks.}
    \item{Type 2:}{
      \eqn{G_1 = g_1 \sqrt{n(n-1)} / (n-2)}{
	G_1 = g_1 * sqrt(n(n-1)) / (n-2)}.
      Used in SAS and SPSS.
    }
    \item{Type 3:}{
      \eqn{b_1 = m_3 / s^3 = g_1 ((n-1)/n)^{3/2}}{
	b_1 = m_3 / s^3 = g_1 ((n-1)/n)^(3/2)}.
      Used in MINITAB and BMDP.}
  }

  All three skewness measures are unbiased under normality.
  
}
\value{
  The estimated skewness of \code{x}.
}
\references{
  D. N. Joanes and C. A. Gill (1998),
  Comparing measures of sample skewness and kurtosis.
  \emph{The Statistician}, \bold{47}, 183--189.
}
\examples{
x <- rnorm(100)
skewness(x)
}
\keyword{univar}
