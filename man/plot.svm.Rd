\name{plot.svm}
\alias{plot.svm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plot svm objects}
\description{
  Generates a scatter plot of the input data of a \code{svm} fit for
  classification models by highlighting the classes and support
  vectors. Optionally, draws a filled contour plot of the class regions.
}
\usage{
plot.svm(x, data, formula, fill = TRUE, grid = 50, slice = list(), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{svm}}
  \item{data}{data to visualize. Should be the same used for fitting.}
  \item{formula}{formula selecting the visualized two dimensions. Only
    needed if more than two input variables are used.}
  \item{fill}{switch indicating whether a contour plot for the class
    regions should be added.}
  \item{grid}{granularity for the contour plot.}
  \item{slice}{a list of named numeric values for the dimensions held
    constant (only needed if more than two variables are
    used). Dimensions not specified are fixed at 0.}
  \item{\dots}{additional graphics parameters passed to
    \code{filled.contour} and \code{plot}.}
}

\author{David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}}
\seealso{\code{\link{svm}}}

\examples{
## a simple example
library(MASS)
data(cats)
m <- svm(Sex~., data = cats)
plot(m, cats)

## more than two variables: fix 2 dimensions
data(iris)
m2 <- svm(Species~., data = iris)
plot(m2, iris, Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))
}
\keyword{neural}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line
\keyword{nonlinear}% __ONLY ONE__ keyword per line
