\name{cmeanscl}
\title{Fuzzy C-Means Clustering}
\usage{
cmeanscl (x, centers, iter.max=100, verbose=FALSE, dist="euclidean",
        method="cmeans", m=2, rate.par = NULL)
}
\alias{cmeanscl}
\arguments{
  \item{x}{Data matrix}
  \item{centers}{Number of clusters or initial values for cluster centers}
  \item{iter.max}{Maximum number of iterations}
  \item{verbose}{If TRUE, make some output during learning}
  \item{dist}{If "euclidean", the mean square error, if "manhattan ",
    the mean absolute error is computed}
  \item{method}{If "cmeans", then we have the cmeans fuzzy clustering method,
    if "ufcl" we have the On-line Update (Unsupervised Fuzzy
    Competitive learning) method} 
  \item{m}{The degree of fuzzification. It is defined for values greater
    than 1}
  \item{rate.par}{The parameter of the learning rate} 
}
\description{
  
  The data given by \code{x} is clustered by the fuzzy kmeans algorithm.
  
  If \code{centers} is a matrix, its rows are taken as the initial cluster
  centers. If \code{centers} is an integer, \code{centers} rows
  of \code{x} are randomly chosen as initial values.
  
  The algorithm stops when the maximum number of iterations (given by
  \code{iter.max}) is reached.

  If \code{verbose} is TRUE, it displays for each iteration the number
  the value of the objective function.

  If \code{dist} is "euclidean", the distance between the cluster center
  and the data points is the Euclidean distance (ordinary kmeans
  algorithm). If "manhattan", the distance between the
  cluster center and the data points is the sum of the absolute values
  of the distances of the coordinates.
  
  If \code{method} is "cmeans", then we have the kmeans fuzzy clustering
  method. If "ufcl" we have the On-line Update (Unsupervised Fuzzy
  Competitive learning) method, which works by performing an update
  directly after each input signal.

  The parameters \code{m} defines the degree of fuzzification. It is
  defined for real values greater than 1 and the bigger it is the more
  fuzzy the membership values of the clustered data points are.
  
  The parameter \code{rate.par} of the learning rate for the "ufcl"
  algorithm which is by default set to rate.par=0.3 and is taking
  real values in (0 , 1).

}
\value{
  \code{cmeanscl} returns an object of class "fclust".
  \item{centers}{The final cluster centers.}
  \item{cluster}{Vector containing the indices of the clusters where
    the data points are assigned to. The maximum membership value of a
    point is considered for partitioning it to a cluster.}
  \item{size}{The number of data points in each cluster.}
  \item{dist}{The distance measure used.}
  \item{m}{The degree of fuzzification.}
  \item{member}{a matrix with the membership values of the data points
    to the clusters.}
  \item{withinss}{Returns the sum of square distances within the
    clusters.} 
  \item{learning}{a list with elements
    \describe{
      \item{ncenters}{The number of the centers,}
      \item{initcenters}{The initial cluster centers,}  
      \item{iter}{The number of iterations performed,} and
      \item{rate.par}{The learning rate for the "ufcl" algorithm.}
  }
}
  \item{call}{Returns a call in which all of the arguments are
    specified by their names.}
 
}
\author{Evgenia Dimitriadou}
\seealso{\code{\link{plot.fclust}}}
\references{
  Nikhil R. Pal, James C. Bezdek, and Richard J. Hathaway.
  Sequential Competitive Learning and the Fuzzy c-Means Clustering
  Algorithms. Neural Networks, Vol. 9, No. 5, pp. 787-796, 1996.
  }
\examples{
# a 2-dimensional example
x<-rbind(matrix(rnorm(100,sd=0.3),ncol=2),
         matrix(rnorm(100,mean=1,sd=0.3),ncol=2))
cl<-cmeanscl(x,2,20,verbose=TRUE,method="cmeans",m=2)
print(cl)
plot(cl,x)   

# a 3-dimensional example
x<-rbind(matrix(rnorm(150,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=1,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=2,sd=0.3),ncol=3))
cl<-cmeanscl(x,6,20,verbose=TRUE,method="cmeans")
plot(cl,x)

}
\keyword{fclust}





