\name{scaclust}
\alias{scaclust}
\title{Fuzzy Clustering using Scatter Matrices}
\usage{
scaclust(x, centers, iter.max=100, verbose=FALSE, method="ad",
         theta = NULL)
}
\arguments{
  \item{x}{The data matrix, where the columns correspond to the
    variables and the rows to the observations.}
  \item{centers}{Number of clusters or initial values for cluster centers}
  \item{iter.max}{Maximum number of iterations}
  \item{verbose}{If \code{TRUE}, make some output during learning}
  \item{method}{If \code{"ad"}, then we have the Adaptive distances method, if
    \code{"mtv"} the Minimum total volume method, if \code{"sand"} the
    Sum of all normalized determinants method and if \code{"mlm"} the
    Maximum likelihood method (Product of Determinants). Abbreviations
    of the method names are also accepted.}
  \item{theta}{A set of constraints for each cluster}
}
\description{
  Four fuzzy clustering methods, namely the Adaptive distances method,
  the Minimum total volume method, the Sum of all normalized
  determinants and the the Maximum likelihood method (Product of
  Determinants) that are based on the calculation of the scatter
  matrices.
}
  
  \details{
  
  The data given by \code{x} is clustered by 4 fuzzy algorithms based on
  the scatter matrices computation.
  
  If \code{centers} is a matrix, its rows are taken as the initial cluster
  centers. If \code{centers} is an integer, \code{centers} rows
  of \code{x} are randomly chosen as initial values.
  
  The algorithm stops when the maximum number of iterations (given by
  \code{iter.max}) is reached.

  If \code{verbose} is \code{TRUE}, it displays for each iteration the number
  the value of the objective function.

  If \code{method} is \code{"ad"}, then we have the Adaptive distances
  method, if \code{"mtv"} the Minimum total volume method, if
  \code{"sand"} the Sum of all normalized determinants method and if
  \code{"mlm"} the Maximum likelihood method (Product of
  Determinants). Note that all these algorithms are adapted for a
  fuzzification parameter of a value \emph{2}.

  \code{theta} is by default \emph{1.0} for every cluster. The relative volumes
  of the clusters are constrained a priori by these constants. An
  inappropriate choice can lead to a bad clustering. The Maximum
  likelihood method does not need this parameter.

}
\value{
  \code{scaclust} returns an object of class \code{"fclust"}.
  \item{centers}{The final cluster centers.}
  \item{size}{The number of data points in each cluster.}
  \item{cluster}{Vector containing the indices of the clusters where
    the data points are assigned to. The maximum membership value of a
    point is considered for partitioning it to a cluster.}
  \item{iter}{The number of iterations performed.}
  \item{membership}{a matrix with the membership values of the data points
    to the clusters.}
  \item{withinerror}{Returns the value of the error function.}
  \item{call}{Returns a call in which all of the arguments are
    specified by their names.}
  
}
\author{Evgenia Dimitriadou}
\references{
  P. J. Rousseeuw, L. Kaufman, and E. Trauwaert. \emph{Fuzzy Clustering using
    Scatter Matrices. Computational Statistics & Data Analysis}, vol.\bold{23},
  p.135-151, 1996.
}
\examples{
## a 2-dimensional example
x<-rbind(matrix(rnorm(100,sd=0.3),ncol=2),
         matrix(rnorm(100,mean=1,sd=0.3),ncol=2))
cl<-scaclust(x,2,20,verbose=TRUE,method="ad")
print(cl)

}
\keyword{cluster}







