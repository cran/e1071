\name{allShortestPaths}
\alias{allShortestPaths}
\alias{extractPath}
\title{Find Shortest Paths Between All Nodes in a Directed Graph}
\description{
  \code{allShortestPaths} finds all shortest paths in a directed (or
  undirected) graph using Floyd's algorithm. \code{extractPath} can be
  used to actually extract the path between a given pair of nodes.
}
\usage{
allShortestPaths(x)
extractPath(obj, start, end)
}
\arguments{
  \item{x}{matrix or distance object}
  \item{obj}{return value of \code{allShortestPaths}}
  \item{start}{integer, starting point of path}
  \item{end}{integer, end point of path}
}
\details{
  If \code{x} is a matrix, then \code{x[i,j]} has to be the length of the
  direct path from point \code{i} to point \code{j}. If no direct
  connection from point \code{i} to point \code{j} exist, then
  \code{x[i,j]} should be either \code{NA} or \code{Inf}. Note that the
  graph can be directed, hence \code{x[i,j]} need not be the same as
  \code{x[j,i]}. The main diagonal of \code{x} is ignored.
  Alternatively, \code{x} can be a distance object as returned by
  \code{\link{dist}} (corresponding to an undirected graph).
}
\value{
  \code{allShortestPaths} returns a list with components
  \item{length}{A matrix with the total lengths of the shortest
    path between each pair of points.}
  \item{middlePoints}{A matrix giving a point in the middle of each
    shortest path (or 0 if the direct connection is the shortest path),
    this is mainly used as input for \code{extractPath}.}
  \code{extractPath} returns a vector of node numbers giving with the
  shortest path between two points.
}
\references{Kumar, V., Grama, A., Gupta, A. and Karypis, G. Introduction
  to Parallel Programming - Design and Analysis of
  Algorithms, Benjamin Cummings Publishing, 1994, ISBN 0-8053-3170-0}
\author{Friedrich Leisch}
\examples{
## build a graph with 5 nodes
x <- matrix(NA, 5, 5)
diag(x) <- 0
x[1,2] <- 30; x[1,3] <- 10
x[2,4] <- 70; x[2,5] <- 40
x[3,4] <- 50; x[3,5] <- 20
x[4,5] <- 60
x[5,4] <- 10
print(x)

## compute all path lengths
z <- allShortestPaths(x)
print(z)

## the following should give 1 -> 3 -> 5 -> 4
extractPath(z, 1, 4)
}
\keyword{optimize}
