\name{plot.fclust}
\title{Plot Clustered Data}
\usage{
plot.fclust(clobj, x, centers=TRUE, initcenters=FALSE)
}
\alias{plot.clobj}
\arguments{
 \item{clobj}{Object returned by a clustering algorithm such as
   \code{\link{cmeanscl}}}
 \item{x}{Data matrix}
 \item{centers}{Mark cluster centers with \code{"X"} (for 2-dimensional
   data only)}
 \item{initcenters}{Mark initial cluster centers with \code{"+"} (for
   2-dimensional data only)}
}
\description{
  Plot the data using different colors for each cluster. If \code{x} has
  more than 2 columns, a \code{\link{pairs}} plot is produced.
}
\author{Friedrich Leisch, Andreas Weingessel and Evgenia Dimitriadou}
\seealso{\code{\link{cmeanscl}},\code{\link{cshell}},\code{\link{scaclust}}}
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
print(cl)
plot(cl,x)
}
\keyword{fclust}







