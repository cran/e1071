\name{ica}
\alias{ica}
\alias{plot.ica}
\alias{print.ica}
\title{Independent Component Analysis}
\usage{
ica(X, lrate, epochs=100, ncomp=dim(X)[2], fun="negative")
}
\arguments{
 \item{X}{The matrix for which the ICA is to be computed}
 \item{lrate}{learning rate}
 \item{epochs}{number of iterations}
 \item{ncomp}{number of independent components}
 \item{fun}{function used for the nonlinear computation part}
}
\description{
This is an R-implementation of the Matlab-Function of
Petteri.Pajunen@hut.fi.

For a data matrix X independent components are extracted by applying a
nonlinear PCA algorithm. The parameter \code{fun} determines which
nonlinearity is used. \code{fun} can either be a function or one of the
following strings "negative kurtosis", "positive kurtosis", "4th
moment" which can be abbreviated to uniqueness. If \code{fun} equals
"negative (positive) kurtosis" the function tanh (x-tanh(x)) is used
which provides ICA for sources with negative (positive) kurtosis. For
\code{fun == "4th moments"} the signed square function is used.
}
\value{
  An object of class \code{"ica"} which is a list with components
  \item{weights}{ICA weight matrix}
  \item{projection}{Projected data}
  \item{epochs}{Number of iterations}
  \item{fun}{Name of the used function}
  \item{lrate}{Learning rate used}
  \item{initweights}{Initial weight matrix}
}

\references{
  Oja et al., ``Learning in Nonlinear Constrained Hebbian Networks'', in
  Proc. ICANN-91, pp. 385--390.

  Karhunen and Joutsensalo, ``Generalizations of Principal Component
  Analysis, Optimization Problems, and Neural Networks'', Neural Networks,
  v. 8, no. 4, pp. 549--562, 1995.
}
\note{Currently, there is no reconstruction from the ICA subspace to the
  original input space.}
\author{Andreas Weingessel}
\keyword{multivariate}
