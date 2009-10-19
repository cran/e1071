\name{matchClasses}
\alias{matchClasses}
\alias{compareMatchedClasses}
\title{Find Similar Classes in Two-way Contingency Tables}
\usage{
matchClasses(tab, method="rowmax", iter=1, maxexact=9, verbose=TRUE)
compareMatchedClasses(x, y, method="rowmax", iter=1,
                      maxexact=9, verbose=FALSE)
}
\arguments{
 \item{tab}{Two-way contingency table of class memberships}
 \item{method}{One of \code{"rowmax"}, \code{"greedy"} or
     \code{"exact"}.}
 \item{iter}{Number of iterations used in greedy search.}
 \item{verbose}{If \code{TRUE}, display some status messages during
   computation.}
 \item{maxexact}{Maximum number of variables for which all possible
   permutations are computed.}
 \item{x, y}{Vectors or matrices with class memberships.}
}
\description{
    Try to find a mapping between the two groupings, such that as many
    cases as possible are in one of the matched pairs. 
}
\details{
    If \code{method="rowmax"}, then each class defining a row in the
    contingency table is mapped to the column of the corresponding row
    maximum. Hence, some columns may be mapped to more than one row
    (while each row is mapped to a single column).

    If \code{method="greedy"} or \code{method="exact"}, then the
    contingency table must be a square matrix and a unique mapping is
    computed. This corresponds to a permutation of columns and rows,
    such that sum of the main diagonal, i.e., the trace of the matrix,
    gets as large as possible. For both methods, first all pairs where
    row and columns maxima correspond and are bigger than the sum of all
    other elements in the corresponding columns and rows together are
    located and fixed (this is a necessary condition for maximal trace).

    If \code{method="exact"}, then for the remaining rows and columns,
    all possible permutations are computed and the optimum is
    returned. This can get computationally infeasible very fast. If more
    than \code{maxexact} rows and columns remain after applying the
    necessary condition, then \code{method} is reset to \code{"greedy"}. If
    \code{method="greedy"}, then a greedy heuristic is tried \code{iter}
    times. Repeatedly a row is picked at random and matched to the free
    column with the maximum value.

    \code{compareMatchedClasses()} computes the contingency table for
    each combination of columns from \code{x} and \code{y} and applies
    \code{matchClasses} to that table. The columns of the table are
    permuted accordingly and then the table is
    passed to \code{\link{classAgreement}}. The resulting agreement
    coefficients (diag, kappa, \ldots) are returned. The return value of
    \code{compareMatchedClasses()} is a list containing a matrix for
    each coefficient; with element (k,l) corresponding to the k-th
    column of \code{x} and l-th column of \code{y}. If \code{y} is
    missing, then the columns of \code{x} are compared with each other.
}
\author{Friedrich Leisch}
\seealso{\code{\link{classAgreement}}}
\examples{
## a stupid example with no class correlations:
g1 <- sample(1:5, size=1000, replace=TRUE)
g2 <- sample(1:5, size=1000, replace=TRUE)
tab <- table(g1, g2)
matchClasses(tab, "exact")

## let pairs (g1=1,g2=4) and (g1=3,g2=1) agree better
k <- sample(1:1000, size=200)
g1[k] <- 1
g2[k] <- 4

k <- sample(1:1000, size=200)
g1[k] <- 3
g2[k] <- 1

tab <- table(g1, g2)
matchClasses(tab, "exact")

## get agreement coefficients:
compareMatchedClasses(g1, g2, method="exact")
}
\keyword{category}
