\name{cshell}
\title{Fuzzy C-Shell Clustering}
\usage{
cshell (x, centers, iter.max=100, verbose=FALSE, dist="euclidean",
        method="cshell", m=2, radius = NULL)

}
\alias{cshell}
\arguments{
  \item{x}{Data matrix}
  \item{centers}{Number of clusters or initial values for cluster centers}
  \item{iter.max}{Maximum number of iterations}
  \item{verbose}{If TRUE, make some output during learning}
  \item{dist}{If "euclidean", the mean square error, if "manhattan ",
    the mean absolute error is computed}
  \item{method}{If "cshell", then we have the c-shell fuzzy clustering method}
  \item{m}{The degree of fuzzification. It is defined for values greater
    than 1}
  \item{radius}{The radius of resulting clusters} 
}
\description{
  
  The data given by \code{x} is clustered by the fuzzy c-shell algorithm.
  
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
  
  If \code{method} is "cshell", then we have the c-shell fuzzy clustering
  method. 

  The parameters \code{m} defines the degree of fuzzification. It is
  defined for real values greater than 1 and the bigger it is the more
  fuzzy the membership values of the clustered data points are.
  
  The parameter \code{radius} is by default set to 0.2 for every cluster.

}
\value{
  \code{cshell} returns an object of class "fclust".
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
      \item{initcenters}{The initial cluster centers,} and 
      \item{iter}{The number of iterations performed.}
      \item{radius}{Returns the radius of the clusters.}
    }
  }
  \item{call}{Returns a call in which all of the arguments are
    specified by their names.}
 
}
\author{Evgenia Dimitriadou}
\seealso{\code{\link{plot.fclust}}}
\references{
  Rajesh N. Dave. Fuzzy Shell-Clustering and Applications to Circle
  Detection in Digital Images. Int. J. of General Systems, Vol. 16,
  pp. 343-355, 1996.
}
\examples{
# a 2-dimensional example
x<-rbind(matrix(rnorm(50,sd=0.3),ncol=2),
         matrix(rnorm(50,mean=1,sd=0.3),ncol=2))
cl<-cshell(x,2,20,verbose=TRUE,method="cshell",m=2)
print(cl)
plot(cl,x)
}
\keyword{fclust}





