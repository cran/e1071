\name{cmeans}
\alias{cmeans}
\alias{print.fclust}
\title{Fuzzy C-Means Clustering}
\usage{
cmeans (x, centers, iter.max=100, verbose=FALSE, dist="euclidean",
        method="cmeans", m=2, rate.par = NULL)
}
\arguments{
  \item{x}{The data matrix where columns correspond to variables and rows to observations}
  \item{centers}{Number of clusters or initial values for cluster centers}
  \item{iter.max}{Maximum number of iterations}
  \item{verbose}{If \code{TRUE}, make some output during learning}
  \item{dist}{Must be one of the following: If \code{"euclidean"}, the
    mean square error, if \code{"manhattan"}, the mean absolute error is
    computed. Abbreviations are also accepted.}
  \item{method}{If \code{"cmeans"}, then we have the cmeans fuzzy
    clustering method, if \code{"ufcl"} we have the On-line Update. Abbreviations in the method names are also accepted.}
  \item{m}{The degree of fuzzification. It is defined for values greater
    than 1}
  \item{rate.par}{The parameter of the learning rate} 
}
\description{
  The fuzzy version of the known \emph{k}means clustering algorithm as
  well as its online update (Unsupervised Fuzzy Competitive learning).  }
\details{
  
  The data given by \code{x} is clustered by the fuzzy \emph{k}means algorithm.
  
  If \code{centers} is a matrix, its rows are taken as the initial cluster
  centers. If \code{centers} is an integer, \code{centers} rows
  of \code{x} are randomly chosen as initial values.
  
  The algorithm stops when the maximum number of iterations (given by
  \code{iter.max}) is reached.

  If \code{verbose} is \code{TRUE}, it displays for each iteration the number
  the value of the objective function.

  If \code{dist} is \code{"euclidean"}, the distance between the
  cluster center and the data points is the Euclidean distance (ordinary
  fuzzy \emph{k}means algorithm). If \code{"manhattan"}, the distance
  between the cluster center and the data points is the sum of the
  absolute values of the distances of the coordinates.
  
  If \code{method} is \code{"cmeans"}, then we have the kmeans fuzzy
  clustering method. If \code{"ufcl"} we have the On-line Update
  (Unsupervised Fuzzy Competitive learning) method, which works by
  performing an update directly after each input signal. 

  The parameters \code{m} defines the degree of fuzzification. It is
  defined for real values greater than 1 and the bigger it is the more
  fuzzy the membership values of the clustered data points are.
  
  The parameter \code{rate.par} of the learning rate for the \code{"ufcl"}
  algorithm which is by default set to \code{rate.par=0.3} and is taking
  real values in (0 , 1).

}
\value{
  \code{cmeans} returns an object of class \code{"fclust"}.
  \item{centers}{The final cluster centers.}
  \item{size}{The number of data points in each cluster.}
  \item{cluster}{Vector containing the indices of the clusters where
    the data points are assigned to. The maximum membership value of a
    point is considered for partitioning it to a cluster.}
  \item{iter}{The number of iterations performed,} and
  \item{membership}{a matrix with the membership values of the data points
    to the clusters.}
  \item{withinerror}{Returns the sum of square distances within the
    clusters.} 
  \item{call}{Returns a call in which all of the arguments are
    specified by their names.}
 
}
\author{Evgenia Dimitriadou}
\references{
  Nikhil R. Pal, James C. Bezdek, and Richard J. Hathaway.
  \emph{Sequential Competitive Learning and the Fuzzy c-Means Clustering
  Algorithms.} Neural Networks, Vol. \bold{9}, No. 5, pp. 787-796, 1996.
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

# assign classes to some new data
y<-rbind(matrix(rnorm(33,sd=0.3),ncol=3),
         matrix(rnorm(33,mean=1,sd=0.3),ncol=3),
         matrix(rnorm(3,mean=2,sd=0.3),ncol=3))
#         ycl<-predict(cl, y, type="both")
}

\keyword{cluster}



