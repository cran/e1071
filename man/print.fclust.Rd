\name{print.fclust}
\alias{print.fclust}
\title{Print Number and Sizes of Clusters}
\usage{
print.fclust(clobj)
}
\arguments{
 \item{clobj}{Object returned by a clustering algorithm such as
   \code{\link{cmeanscl}}.}
}
\description{
  Prints the number of the clusters and the sizes of them.
}
\author{Friedrich Leisch, Andreas Weingessel and Evgenia Dimitriadou}
\seealso{
  \code{\link{cmeanscl}},
  \code{\link{cshell}},
  \code{\link{scaclust}}
}
\examples{
## a 2-dimensional example
x<-rbind(matrix(rnorm(100,sd=0.3),ncol=2),
         matrix(rnorm(100,mean=1,sd=0.3),ncol=2))
cl<-cmeanscl(x,2,20,verbose=TRUE,method="cmeans",m=2)
print(cl)

## a 3-dimensional example
x<-rbind(matrix(rnorm(150,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=1,sd=0.3),ncol=3),
         matrix(rnorm(150,mean=2,sd=0.3),ncol=3))
cl<-cmeanscl(x,6,20,verbose=TRUE,method="cmeans")
print(cl)
}
\keyword{cluster}
