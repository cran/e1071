\name{bootstrap.lca}

\alias{bootstrap.lca}
\alias{print.bootstrap.lca}

\title{Bootstrap Samples of LCA Results}
\description{
This function draws bootstrap samples from a given LCA model and refits
a new LCA model for each sample. The quality of fit of these models is
compared to the original model.
}
\usage{
bootstrap.lca(l, nsamples=10, lcaiter=30, verbose=FALSE)
}

\arguments{
  \item{l}{An LCA model as created by \code{\link{lca}}}
  \item{nsamples}{Number of bootstrap samples}
  \item{lcaiter}{Number of LCA iterations}
  \item{verbose}{If \code{TRUE} some output is printed during the
    computations.} 
}
\details{
From a given LCA model \code{l}, \code{nsamples} bootstrap samples are
drawn. For each sample a new LCA model is fitted. The goodness of fit
for each model is computed via Likelihood Ratio and Pearson's
Chisquare. The values for the fitted models are compared with the values
of the original model \code{l}. By this method it can be tested whether
the data to which \code{l} was originally fitted come from an LCA model.
}
\value{
  An object of class \code{bootstrap.lca} is returned, containing
  \item{logl, loglsat}{The LogLikelihood of the models and of the
    corresponding saturated models}
  \item{lratio}{Likelihood quotient of the models and the corresponding
    saturated models}
  \item{lratiomean, lratiosd}{Mean and Standard deviation of
    \code{lratio}}
  \item{lratioorg}{Likelihood quotient of the original model and the
    corresponding saturated model}
  \item{zratio}{Z-Statistics of \code{lratioorg}}
  \item{pvalzratio, pvalratio}{P-Values for \code{zratio}, computed via normal
    distribution and empirical distribution}
  \item{chisq}{Pearson's Chisq of the models}
  \item{chisqmean, chisqsd}{Mean and Standard deviation of
    \code{chisq}}
  \item{chisqorg}{Pearson's Chisq of the original model}
  \item{zchisq}{Z-Statistics of \code{chisqorg}}
  \item{pvalzchisq, pvalchisq}{P-Values for \code{zchisq}, computed via normal
    distribution and empirical distribution}
  \item{nsamples}{Number of bootstrap samples}
  \item{lcaiter}{Number of LCA Iterations}
}
\references{Anton K. Formann: ``Die Latent-Class-Analysis'', Beltz
    Verlag 1984}

\author{Andreas Weingessel}
    
\seealso{\code{\link{lca}}}

\examples{
## Generate a 4-dim. sample with 2 latent classes of 500 data points each.
## The probabilities for the 2 classes are given by type1 and type2.
type1 <- c(0.8,0.8,0.2,0.2)
type2 <- c(0.2,0.2,0.8,0.8)
x <- matrix(runif(4000),nr=1000)
x[1:500,] <- t(t(x[1:500,])<type1)*1
x[501:1000,] <- t(t(x[501:1000,])<type2)*1

l <- lca(x, 2, niter=5)
bl <- bootstrap.lca(l,nsamples=3,lcaiter=5)
bl
}

\keyword{multivariate}

