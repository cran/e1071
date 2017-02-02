\name{impute}
\alias{impute}
\title{Replace Missing Values}
\usage{
impute(x, what = c("median", "mean"))
}
\arguments{
 \item{x}{A matrix or dataframe.}
 \item{what}{What to impute.}
}
\description{
    Replaces missing values of a matrix or dataframe with the
    medians (\code{what="median"}) or means (\code{what="mean"}) of the
    respective columns.}
\value{
    A matrix or dataframe.
}
\author{Friedrich Leisch}
\examples{
x<- matrix(1:10, ncol=2)
x[c(1,3,7)] <- NA
print(x)
print(impute(x))
}
\keyword{manip}
