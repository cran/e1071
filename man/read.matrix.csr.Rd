\name{read.matrix.csr}
\alias{read.matrix.csr}
\alias{write.matrix.csr}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Read/Write Sparse Data}
\description{
  reads and writes a file in sparse data format.
}
\usage{
read.matrix.csr(file, fac = TRUE, ncol = NULL)
write.matrix.csr(x, file = "out.dat", y = NULL, fac = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{matrix.csr}}
  \item{y}{A vector (either numeric or a factor)}
  \item{file}{The filename.}
  \item{fac}{If \code{TRUE}, the
    y-values (if any) are interpreted as factor levels.}
  \item{ncol}{Number of columns, detected automatically. Can be used to
    add empty columns (possibly not stored in the sparse format).}
}
\value{
  If the data file includes no y variable, \code{read.matrix.csr} returns an object of class
  \code{matrix.csr}, else a list with components:
  \item{x}{object of class \code{matrix.csr}}
  \item{y}{vector of numeric values or factor levels, depending on \code{fac}.}
}
\author{
  David Meyer\cr
  \email{David.Meyer@R-project.org}
}
\seealso{\code{\link[SparseM]{matrix.csr}}}

\examples{
\dontrun{
library(methods)
if (require(SparseM)) {
    data(iris)
    x <- as.matrix(iris[,1:4])
    y <- iris[,5]
    xs <- as.matrix.csr(x)
    write.matrix.csr(xs, y = y, file = "iris.dat")
    xs2 <- read.matrix.csr("iris.dat")$x
    if (!all(as.matrix(xs) == as.matrix(xs2)))
        stop("Error: objects are not equal!")
}
}
}
\keyword{IO}% at least one, from doc/KEYWORDS
