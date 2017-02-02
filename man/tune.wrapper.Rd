\name{tune.wrapper}
\alias{tune.wrapper}
\alias{tune.rpart}
\alias{best.rpart}
\alias{tune.svm}
\alias{best.svm}
\alias{tune.nnet}
\alias{best.nnet}
\alias{tune.randomForest}
\alias{best.randomForest}
\alias{tune.knn}
\title{Convenience Tuning Wrapper Functions}
\description{
  Convenience tuning wrapper functions, using \code{tune}.
}
\usage{
tune.svm(x, y = NULL, data = NULL, degree = NULL, gamma = NULL, coef0 = NULL,
         cost = NULL, nu = NULL, class.weights = NULL, epsilon = NULL, ...)
best.svm(x, tunecontrol = tune.control(), ...)
 
tune.nnet(x, y = NULL, data = NULL, size = NULL, decay = NULL,
          trace = FALSE, tunecontrol = tune.control(nrepeat = 5), 
          ...)
best.nnet(x, tunecontrol = tune.control(nrepeat = 5), ...)

tune.rpart(formula, data, na.action = na.omit, minsplit = NULL,
           minbucket = NULL, cp = NULL, maxcompete = NULL, maxsurrogate = NULL,
           usesurrogate = NULL, xval = NULL, surrogatestyle = NULL, maxdepth =
           NULL, predict.func = NULL, ...)
best.rpart(formula, tunecontrol = tune.control(), ...)

tune.randomForest(x, y = NULL, data = NULL, nodesize = NULL, 
                  mtry = NULL, ntree = NULL, ...)
best.randomForest(x, tunecontrol = tune.control(), ...)

tune.knn(x, y, k = NULL, l = NULL, ...) 

}
\arguments{
  \item{formula, x, y, data}{formula and data arguments of function to be tuned.}
  \item{predict.func}{predicting function.}
  \item{na.action}{function handling missingness.}
  \item{minsplit, minbucket, cp, maxcompete,
        maxsurrogate, usesurrogate, xval,
        surrogatestyle, maxdepth}{\code{rpart} parameters.}
  \item{degree, gamma, coef0, cost, nu, class.weights, epsilon}{\code{svm}
    parameters.}
  \item{k, l}{\code{knn} parameters.}
  \item{mtry, nodesize, ntree}{\code{randomForest} parameters.}
  \item{size, decay, trace}{parameters passed to
    \code{nnet}.}
  \item{tunecontrol}{object of class \code{"tune.control"} containing
    tuning parameters.}
  \item{\dots}{Further parameters passed to \code{tune}.}
}

\value{
  \code{tune.foo()} returns a tuning object including the best parameter set obtained
  by optimizing over the specified parameter vectors. \code{best.foo()}
  directly returns the best model, i.e. the fit of a new model using the
  optimal parameters found by \code{tune.foo}.
}
\details{For examples, see the help page of \code{tune()}.}

\author{
  David Meyer\cr
  \email{David.Meyer@R-project.org}
}
\seealso{\code{\link{tune}}}

\keyword{models}
