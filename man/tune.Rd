\name{tune}
\alias{tune}
\alias{print.tune}
\alias{summary.tune}
\alias{print.summary.tune}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Parameter tuning of fuctions using grid search}
\description{
  This generic function tunes hyperparameters of statistical methods
  using a grid search over supplied parameter ranges.
}
\usage{
tune(method, train.x, train.y = NULL, data = list(), validation.x =
     NULL, validation.y = NULL, ranges, random = FALSE, nrepeat = 1,
     repeat.aggregate = min, sampling = c("cross", "fix", "bootstrap"),
     sampling.aggregate = mean, cross = 10, fix = 2/3, nboot = 10,
     boot.size = 9/10, predict.func = predict, best.model = TRUE,
     performances = TRUE, ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{method}{function to be tuned.}
  \item{train.x}{either a formula or a matrix of predictors.}
  \item{train.y}{the response variable if \code{train.x} is a predictor
    matrix. Ignored if \code{train.x} is a formula.}
  \item{data}{data, if a formula interface is used. Ignored, if
    predictor matrix and response are supplied directly.}
  \item{validation.x}{an optional validation set. Depending on whether a
    formula interface is used or not, the response can be
    included in \code{validation.x} or separately speciefied using \code{validation.y}.}
  \item{validation.y}{if no formula interface is used, the response of
    the (optional) validation set.}
  \item{ranges}{a named list of parameter vectors spanning the sampling
    space. The vectors will usually be created by \code{seq}.}
  \item{random}{if an integer value is specified, \code{random}
    parameter vectors are drawn from the parameter space.}
  \item{nrepeat}{specifies how often training shall be repeated.}
  \item{repeat.aggregate}{function for aggregating the repeated training results.}
  \item{sampling}{sampling scheme. If \code{sampling = "cross"}, a
    \code{cross}-times cross validation is performed. If \code{sampling
      = "boot"}, \code{nboot} training sets of size \code{boot.size} (part)
    are sampled from the supplied data. If \code{sampling
      = "fix"}, a single split into training/validation set is
    used, the training set containing a \code{fix} part of the supplied
    data. Note that a separate validation set can be supplied via
    \code{validation.x} and \code{validation.y}. It is only used for
    \code{sampling = "boot"} and \code{sampling = "fix"}; in the latter
    case, \code{fix} is set to 1.}
  \item{sampling.aggregate}{function for aggregating the training
    results on the generated training samples.}
  \item{cross}{number of partitions for cross-validation.}
  \item{fix}{part of the data used for training in fixed sampling.}
  \item{nboot}{number of bootstrap replications.}
  \item{boot.size}{size of the bootstrap samples.}
  \item{predict.func}{optional predict function, if the standard \code{predict}
  behaviour is inadequate.}
  \item{best.model}{if \code{TRUE}, the best model is trained and
    returned (the best parameter set is used for
    training on the complete training set).}
  \item{performances}{if \code{TRUE}, the performance results for all
    parameter combinations are returned.}
  \item{\dots}{Further parameters passed to the training functions.}
}
\value{
  An object of class \code{tune}, including the components:
  \item{best.parameters}{a 1 x k data frame, k number of parameters.}
  \item{best.performance}{best achieved performance.}
  \item{performances}{if requested, a data frame of all parameter
    combinations along with the corresponding performance results.}
  \item{if requested, the model trained on the complete training data
    using the best parameter combination.}
}
\note{
  As performance measure, the classification error is used
  for classification, and the mean squared error for regression.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{\code{\link{plot.tune}}, \code{\link{tune.svm}}}
\examples{
  data(iris)
  ## tune `svm' for classification with RBF-kernel (default in svm),
  ## using one split for training/validation set
  
  obj <- tune(svm, Species~., data = iris, sampling = "fix",
              ranges = list(gamma = 2^(-1:1), cost = 2^(2:4))
             )

  ## alternatively:
  ## obj <- tune.svm(Species~., data = iris, gamma = 2^(-1:1), cost = 2^(2:4))

  summary(obj)
  plot(obj)

  ## tune `knn' using a convenience function; this time with the
  ## conventional interface and bootstrap sampling:
  x <- iris[,-5]
  y <- iris[,5]
  obj2 <- tune.knn(x, y, k = 1:5, sampling = "boot")
  summary(obj2)
  plot(obj2)

  ## tune `rpart' for regression, using 10-fold cross validation (default)
  data(mtcars)
  obj3 <- tune.rpart(mpg~., data = mtcars, minsplit = c(5,10,15))
  summary(obj3)
  plot(obj3)
}
\keyword{models}% at least one, from doc/KEYWORDS

