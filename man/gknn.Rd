\name{gknn}
\alias{gknn}
\alias{gknn.default}
\alias{gknn.formula}
\alias{print.gknn}
\alias{predict.gknn}
\title{Generalized k-Nearest Neighbors Classification or Regression}
\description{
\code{gknn} is an implementation of the k-nearest neighbours algorithm making use of general distance measures. A formula interface is provided.
}
\usage{
\method{gknn}{formula}(formula, data = NULL, ..., subset, na.action = na.pass, scale = TRUE)
\method{gknn}{default}(x, y, k = 1, method = NULL, 
                       scale = TRUE, use_all = TRUE, 
                       FUN = mean, ...)
\method{predict}{gknn}(object, newdata, 
                         type = c("class", "votes", "prob"), 
                         ...,
                         na.action = na.pass)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          \sQuote{gknn} is called from.}
  \item{x}{a data matrix.}
  \item{y}{a response vector with one label for each row/component of
    \code{x}. Can be either a factor (for classification tasks)
    or a numeric vector (for regression).}
  \item{k}{number of neighbours considered.}
  \item{scale}{a logical vector indicating the variables to be
    scaled. If \code{scale} is of length 1, the value is recycled as
    many times as needed.
    By default, numeric \emph{matrices} are scaled to zero mean and unit variance. The center and scale
    values are returned and used for later predictions.
    Note that the default metric for data frames is the Gower metric
    which \emph{standardizes} the values to the unit interval.}
  \item{method}{Argument passed to \code{dist()} from the \code{proxy} package to select the distance metric used: a function, or a mnemonic string referencing the distance measure. Defaults to \code{"Euclidean"} for metric matrices, to \code{"Jaccard"} for logical matrices and to \code{"Gower"} for data frames.
  }
  \item{use_all}{controls handling of ties. If true, all distances equal to the kth largest are included. If false, a random selection of distances equal to the kth is chosen to use exactly k neighbours.}
  \item{FUN}{function used to aggregate the k nearest target values in case of regression.}
  \item{object}{object of class \code{gknn}.}
  \item{newdata}{matrix or data frame with new instances.}
  \item{type}{character specifying the return type in case of class
  predictions: for \code{"class"}, the class labels; for \code{"prob"}, the class distribution for all k neighbours considered; for \code{"votes"}, the raw counts.}
  \item{\dots}{additional parameters passed to \code{dist()}}
  \item{subset}{An index vector specifying the cases to be used in the
          training sample.  (NOTE: If given, this argument must be
          named.)}
  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.pass}. (NOTE: If given, this argument must be named.)}	
}

\value{
  For \code{gknn()}, an object of class \code{"gknn"} containing the data and the specified parameters. For \code{predict.gknn()}, a vector of predictions, or a matrix with votes for all classes. In case of an overall class tie, the predicted class is chosen by random.
}
\author{
  David Meyer (\email{David.Meyer@R-project.org})
}
\seealso{
  \code{\link[proxy]{dist}} (in package \pkg{proxy})
}
\examples{
data(iris)

model <- gknn(Species ~ ., data = iris)
predict(model, iris[c(1, 51, 101),])

test = c(45:50, 95:100, 145:150)

model <- gknn(Species ~ ., data = iris[-test,], k = 3, method = "Manhattan")
predict(model, iris[test,], type = "votes")

model <- gknn(Species ~ ., data = iris[-test], k = 3, method = "Manhattan")
predict(model, iris[test,], type = "prob")

}

\keyword{nonlinear}
\keyword{classif}









