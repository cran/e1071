\name{boxplot.bclust}
\alias{boxplot.bclust}
\title{Boxplot of Cluster Profiles}
\usage{
\method{boxplot}{bclust}(x, n=nrow(x$centers), bycluster=TRUE,
        main=deparse(substitute(x)), oneplot=TRUE,
        which=1:n, ...)
}
\arguments{
 \item{x}{Clustering result, object of class \code{"bclust"}.}% 
 \item{n}{Number of clusters to plot, by default the number of clusters
     used in the call of \code{\link{bclust}}.}
 \item{bycluster}{If \code{TRUE} (default), a boxplot for each cluster
     is plotted. If \code{FALSE}, a boxplot for each variable is plotted.}
 \item{main}{Main title of the plot, by default the name of the cluster
     object.} 
 \item{oneplot}{If \code{TRUE}, all boxplots appear on one screen (using
   an appropriate rectangular layout).}
 \item{which}{Number of clusters which should be plotted, default is all
   clusters.}
 \item{...}{Additional arguments for \code{\link{boxplot}}.}
}
\description{
    Makes boxplots of the results of a bagged clustering run.
}
\author{Friedrich Leisch}
\keyword{hplot}
\examples{
data(iris)
bc1 <- bclust(iris[,1:4], 3, base.centers=5)
boxplot(bc1)
}

