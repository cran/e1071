\name{tune.rpart}
\alias{tune.rpart}
\alias{best.rpart}
\alias{tune.svm}
\alias{best.svm}
\alias{tune.nnet}
\alias{best.nnet}
\alias{tune.randomForest}
\alias{best.randomForest}
\alias{rpart.wrapper}
\alias{tune.knn}
\alias{knn.wrapper}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convenience tuning functions}
\description{
  Convenience tuning functions, using \code{tune}.
}
\usage{
tune.svm(x, y = NULL, degree = NULL, gamma = NULL, coef0 = NULL, cost = NULL, 
         nu = NULL, ...)
best.svm(x, ...)
 
tune.nnet(x, y = NULL, size = NULL, decay = NULL, nrepeat = 5, trace = FALSE, 
          predict.func = function(...) predict(..., type = "class"), ...)
best.nnet(x, ...)

tune.rpart(formula, data, na.action = na.omit, minsplit = NULL,
           minbucket = NULL, cp = NULL, maxcompete = NULL, maxsurrogate = NULL,
           usesurrogate = NULL, xval = NULL, surrogatestyle = NULL, maxdepth =
           NULL, predict.func = NULL, ...)
best.rpart(formula, ...)
rpart.wrapper(formula, minsplit=20, minbucket=round(minsplit/3), cp=0.01, 
              maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10,
              surrogatestyle=0, maxdepth=30, ...)

tune.randomForest(x, y = NULL, nodesize = NULL, mtry = NULL, ntree = NULL, ...)
best.randomForest(x, ...)

tune.knn(x, y, k = NULL, l = NULL, ...) 
knn.wrapper(x, y, k = 1, l = 0, ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula, x, y, data}{formula and data arguments of function to be tuned.}
  \item{predict.func}{predicting function.}
  \item{na.action}{function handling missingness.}
  \item{minsplit, minbucket, cp, maxcompete,
        maxsurrogate, usesurrogate, xval,
        surrogatestyle, maxdepth}{\code{rpart} parameters.}
  \item{degree, gamma, coef0, cost, nu}{\code{svm}
    parameters.}
  \item{k, l}{\code{knn} parameters.}
  \item{mtry, nodesize, ntree}{\code{randomForest} parameters.}
  \item{size, decay, nrepeat, trace}{\code{nnet} parameters.}
  \item{\dots}{Further parameters passed to \code{tune}.}
}

\value{
  \code{tune.foo} returns a tuning object including the best parameter set obtained
  by optimizing over the specified parameter vectors. \code{best.foo}
  directly returns the best model, i.e. the fit of a new model using the
  optimal parameters found by \code{tune.foo}.
}

\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{\code{\link{tune}}}

\keyword{models}% __ONLY ONE__ keyword per line
