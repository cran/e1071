\name{cshell}
\alias{cshell}
\title{Fuzzy C-Shell Clustering}
\usage{
cshell(x, centers, iter.max=100, verbose=FALSE, dist="euclidean",
       method="cshell", m=2, radius = NULL)
}
\arguments{
  \item{x}{The data matrix, were columns correspond to the variables and
    rows to observations.}
  \item{centers}{Number of clusters or initial values for cluster centers}
  \item{iter.max}{Maximum number of iterations}
  \item{verbose}{If \code{TRUE}, make some output during learning}
  \item{dist}{Must be one of the following: If \code{"euclidean"}, the
    mean square error, if \code{"manhattan"}, the mean absolute error is
    computed. Abbreviations are also accepted.}
  \item{method}{Currently, only the \code{"cshell"} method; the c-shell fuzzy
    clustering method}
  \item{m}{The degree of fuzzification. It is defined for values greater
    than \emph{1}}
  \item{radius}{The radius of resulting clusters} 
}
\description{
  The \emph{c}-shell clustering algorithm, the shell prototype-based version
  (ring prototypes) of the fuzzy \emph{k}means clustering method.
}
\details{
  
  The data given by \code{x} is clustered by the fuzzy \emph{c}-shell algorithm.
  
  If \code{centers} is a matrix, its rows are taken as the initial cluster
  centers. If \code{centers} is an integer, \code{centers} rows
  of \code{x} are randomly chosen as initial values.
  
  The algorithm stops when the maximum number of iterations (given by
  \code{iter.max}) is reached.

  If \code{verbose} is \code{TRUE}, it displays for each iteration the number
  the value of the objective function.

  If \code{dist} is \code{"euclidean"}, the distance between the
  cluster center and the data points is the Euclidean distance (ordinary
  kmeans algorithm). If \code{"manhattan"}, the distance between the
  cluster center and the data points is the sum of the absolute values
  of the distances of the coordinates.
  
  If \code{method} is \code{"cshell"}, then we have the \emph{c}-shell
  fuzzy clustering method.

  The parameters \code{m} defines the degree of fuzzification. It is
  defined for real values greater than 1 and the bigger it is the more
  fuzzy the membership values of the clustered data points are.
  
  The parameter \code{radius} is by default set to \emph{0.2} for every
  cluster.

}
\value{
  \code{cshell} returns an object of class \code{"cshell"}.
  \item{centers}{The final cluster centers.}
  \item{size}{The number of data points in each cluster.}
  \item{cluster}{Vector containing the indices of the clusters where
    the data points are assigned to. The maximum membership value of a
    point is considered for partitioning it to a cluster.}
  \item{iter}{The number of iterations performed.}
  \item{membership}{a matrix with the membership values of the data points
    to the clusters.}
  \item{withinerror}{Returns the sum of square distances within the
    clusters.} 
  \item{call}{Returns a call in which all of the arguments are
    specified by their names.}
  
}
\author{Evgenia Dimitriadou}
\references{
  Rajesh N. Dave. \emph{Fuzzy Shell-Clustering and Applications to Circle
  Detection in Digital Images.} Int. J. of General Systems, Vol. \bold{16},
  pp. 343-355, 1996.
}
\examples{
## a 2-dimensional example
x<-rbind(matrix(rnorm(50,sd=0.3),ncol=2),
         matrix(rnorm(50,mean=1,sd=0.3),ncol=2))
cl<-cshell(x,2,20,verbose=TRUE,method="cshell",m=2)
print(cl)

# assign classes to some new data
y<-rbind(matrix(rnorm(13,sd=0.3),ncol=2),
         matrix(rnorm(13,mean=1,sd=0.3),ncol=2))
#         ycl<-predict(cl, y, type="both")
        }
\keyword{cluster}
