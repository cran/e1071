\name{lca}

\alias{lca}
\alias{print.lca}
\alias{summary.lca}
\alias{print.summary.lca}
\alias{predict.lca}

\title{Latent Class Analysis (LCA)}
\usage{
l <- lca(x, nclass, niter=100, matchdata=FALSE, verbose=FALSE)
print(l)
summary(l)
predict(l,x)
}

\arguments{
 \item{x}{Either a data matrix of binary observations or a list of
     patterns as created by \code{\link{countpattern}}}
 \item{nclass}{Number of classes used for LCA}
 \item{niter}{Number of Iterations}
 \item{matchdata}{If TRUE and \code{x} is a data matrix, the class
     membership of every data point is returned, otherwise the class
     membership of every pattern is returned.}
 \item{verbose}{If TRUE some output is printed during the computations.}
}

\description{
A latent class analysis with \code{nclass} classes is performed on the
data given by \code{x}.
}

\value{
    An object of class \code{lca} is returned, containing
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

\seealso{\code{\link{countpattern}}, \code{\link{bootstrap.lca}}}

%\examples{
%## No examples with small data sets

%library(sfb)
%data(scen.1)
%l <- lca(scen.1$x, 6, niter=20)
%print(l)
%summary(l)
%p <- predict(l, scen.1$x)
%table(p, scen.1$classes)
%}

\keyword{multivariate, cluster}
