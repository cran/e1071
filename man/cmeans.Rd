\name{cmeans}
\alias{cmeans}
\alias{print.fclust}
\title{Fuzzy C-Means Clustering}
\description{
  The fuzzy version of the known \emph{k}means clustering algorithm as
  well as an on-line variant (Unsupervised Fuzzy Competitive learning).
}
\usage{
cmeans(x, centers, iter.max = 100, verbose = FALSE,
       dist = "euclidean", method = "cmeans", m = 2,
       rate.par = NULL, weights = 1, control = list())
}
\arguments{
  \item{x}{The data matrix where columns correspond to variables and
    rows to observations.} 
  \item{centers}{Number of clusters or initial values for cluster
    centers.}
  \item{iter.max}{Maximum number of iterations.}
  \item{verbose}{If \code{TRUE}, make some output during learning.}
  \item{dist}{Must be one of the following: If \code{"euclidean"}, the
    mean square error, if \code{"manhattan"}, the mean absolute error is
    computed.  Abbreviations are also accepted.}
  \item{method}{If \code{"cmeans"}, then we have the \eqn{c}-means fuzzy
    clustering method, if \code{"ufcl"} we have the on-line update.
    Abbreviations are also accepted.}
  \item{m}{A number greater than 1 giving the degree of fuzzification.}
  \item{rate.par}{A number between 0 and 1 giving the parameter of the
    learning rate for the on-line variant.  The default corresponds to
    \eqn{0.3}.}
  \item{weights}{a numeric vector with non-negative case weights.
    Recycled to the number of observations in \code{x} if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}

}
\details{
  The data given by \code{x} is clustered by generalized versions of the
  fuzzy \emph{c}-means algorithm, which use either a fixed-point or an
  on-line heuristic for minimizing the objective function
  \deqn{\sum_i \sum_j w_i u_{ij}^m d_{ij},}
  where \eqn{w_i} is the weight of observation \eqn{i}, \eqn{u_{ij}} is
  the membership of observation \eqn{i} in cluster \eqn{j}, and
  \eqn{d_{ij}} is the distance (dissimilarity) between observation
  \eqn{i} and center \eqn{j}.  The dissimilarities used are the sums of
  squares (\code{"euclidean"}) or absolute values (\code{"manhattan"})
  of the element-wise differences.
  
  If \code{centers} is a matrix, its rows are taken as the initial cluster
  centers.  If \code{centers} is an integer, \code{centers} rows of
  \code{x} are randomly chosen as initial values.
  
  The algorithm stops when the maximum number of iterations (given by
  \code{iter.max}) is reached, or when the algorithm is unable to reduce
  the current value \code{val} of the objective function by
  \code{reltol * (abs(val) * reltol)} at a step.  The relative
  convergence tolerance \code{reltol} can be specified as the
  \code{reltol} component of the list of control parameters, and
  defaults to \code{sqrt(.Machine$double.eps)}.

  If \code{verbose} is \code{TRUE}, each iteration displays its number
  and the value of the objective function.

  If \code{method} is \code{"cmeans"}, then we have the \eqn{c}-means
  fuzzy clustering method, see for example Bezdek (1981).  If
  \code{"ufcl"}, we have the On-line Update (Unsupervised Fuzzy
  Competitive Learning) method due to Chung and Lee (1992), see also Pal
  et al (1996).  This method works by performing an update directly
  after each input signal (i.e., for each single observation).

  The parameters \code{m} defines the degree of fuzzification.  It is
  defined for real values greater than 1 and the bigger it is the more
  fuzzy the membership values of the clustered data points are.
}
\value{
  An object of class \code{"fclust"} which is a list with components:
  \item{centers}{the final cluster centers.}
  \item{size}{the number of data points in each cluster of the closest
    hard clustering.}
  \item{cluster}{a vector of integers containing the indices of the
    clusters where the data points are assigned to for the closest hard
    clustering, as obtained by assigning points to the (first) class with
    maximal membership.}
  \item{iter}{the number of iterations performed.}
  \item{membership}{a matrix with the membership values of the data points
    to the clusters.}
  \item{withinerror}{the value of the objective function.}
  \item{call}{the call used to create the object.}
}
\author{
  Evgenia Dimitriadou and Kurt Hornik
}
\references{
  J. C. Bezdek (1981).
  \emph{Pattern recognition with fuzzy objective function algorithms}.
  New York: Plenum.

  Fu Lai Chung and Tong Lee (1992).
  Fuzzy competitive learning.
  \emph{Neural Networks}, \bold{7}(3), 539--551.

  Nikhil R. Pal, James C. Bezdek, and Richard J. Hathaway (1996).
  Sequential competitive learning and the fuzzy c-means clustering
  algorithms.
  \emph{Neural Networks}, \bold{9}(5), 787--796.
}
\examples{
# a 2-dimensional example
x<-rbind(matrix(rnorm(100,sd=0.3),ncol=2),
         matrix(rnorm(100,mean=1,sd=0.3),ncol=2))
cl<-cmeans(x,2,20,verbose=TRUE,method="cmeans",m=2)
print(cl)

# a 3-dimensional example
x<-rbind(matrix(rnorm(150,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=1,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=2,sd=0.3),ncol=3))
cl<-cmeans(x,6,20,verbose=TRUE,method="cmeans")
print(cl)
}
\keyword{cluster}
