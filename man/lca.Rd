\name{lca}
\alias{lca}
\alias{print.lca}
\alias{summary.lca}
\alias{print.summary.lca}
\alias{predict.lca}
\title{Latent Class Analysis (LCA)}
\usage{
lca(x, k, niter=100, matchdata=FALSE, verbose=FALSE)
}
\arguments{
 \item{x}{Either a data matrix of binary observations or a list of
     patterns as created by \code{\link{countpattern}}}
 \item{k}{Number of classes used for LCA}
 \item{niter}{Number of Iterations}
 \item{matchdata}{If \code{TRUE} and \code{x} is a data matrix, the class
     membership of every data point is returned, otherwise the class
     membership of every pattern is returned.}
 \item{verbose}{If \code{TRUE} some output is printed during the
   computations.}
}
\description{
  A latent class analysis with \code{k} classes is performed on the data
  given by \code{x}.
}
\value{
  An object of class \code{"lca"} is returned, containing
  \item{w}{Probabilities to belong to each class}
  \item{p}{Probabilities of a `1' for each variable in each class}
  \item{matching}{Depending on \code{matchdata} either the class
      membership of each pattern or of each data point}
  \item{logl, loglsat}{The LogLikelihood of the model and of the
      saturated model}
  \item{bic, bicsat}{The BIC of the model and of the
      saturated model}
  \item{chisq}{Pearson's Chisq}
  \item{lhquot}{Likelihood quotient of the model and the saturated
      model}
  \item{n}{Number of data points.}
  \item{np}{Number of free parameters.}
}

\references{Anton K. Formann: ``Die Latent-Class-Analysis'', Beltz
    Verlag 1984}

\author{Andreas Weingessel}

\seealso{
  \code{\link{countpattern}},
  \code{\link{bootstrap.lca}}
}
\examples{
## Generate a 4-dim. sample with 2 latent classes of 500 data points each.
## The probabilities for the 2 classes are given by type1 and type2.
type1 <- c(0.8,0.8,0.2,0.2)
type2 <- c(0.2,0.2,0.8,0.8)
x <- matrix(runif(4000),nr=1000)
x[1:500,] <- t(t(x[1:500,])<type1)*1
x[501:1000,] <- t(t(x[501:1000,])<type2)*1

l <- lca(x, 2, niter=5)
print(l)
summary(l)
p <- predict(l, x)
table(p, c(rep(1,500),rep(2,500)))
}
\keyword{multivariate}
\keyword{cluster}
