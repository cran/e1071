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
     NULL, validation.y = NULL, ranges = NULL, predict.func = predict,
     control = tune.control(), ...)
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
  \item{predict.func}{optional predict function, if the standard \code{predict}
    behaviour is inadequate.}
  \item{control}{object of class \code{"tune.control"}, as created by the
    function \code{tune.control()}.}
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
\details{
  As performance measure, the classification error is used
  for classification, and the mean squared error for regression. It is
  possible to specify only one parameter combination (i.e., vectors of
  length 1) to obtain an error estimation of the specified type
  (bootstrap, cross-classification, etc.) on the given data set. For
  conveneince, there
  are several \code{tune.foo()} wrappers defined, e.g., for
  \code{nnet()}, \code{randomForest()},
  \code{rpart()}, \code{svm()}, and \code{knn()}.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{\code{\link{tune.control}}, \code{\link{plot.tune}},
  \code{\link{tune.svm}}, \link{tune.wrapper}}
\examples{
  data(iris)
  ## tune `svm' for classification with RBF-kernel (default in svm),
  ## using one split for training/validation set
  
  obj <- tune(svm, Species~., data = iris, 
              ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
              control = tune.control(sampling = "fix")
             )

  ## alternatively:
  ## obj <- tune.svm(Species~., data = iris, gamma = 2^(-1:1), cost = 2^(2:4))

  summary(obj)
  plot(obj)

  ## tune `knn' using a convenience function; this time with the
  ## conventional interface and bootstrap sampling:
  x <- iris[,-5]
  y <- iris[,5]
  obj2 <- tune.knn(x, y, k = 1:5, control = tune.control(sampling = "boot"))
  summary(obj2)
  plot(obj2)

  ## tune `rpart' for regression, using 10-fold cross validation (default)
  data(mtcars)
  obj3 <- tune.rpart(mpg~., data = mtcars, minsplit = c(5,10,15))
  summary(obj3)
  plot(obj3)

  ## simple error estimation for lm using 10-fold cross validation
  tune(lm, mpg~., data = mtcars)
}
\keyword{models}% at least one, from doc/KEYWORDS

