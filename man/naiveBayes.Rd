\name{naiveBayes}
\alias{naiveBayes}
\alias{naiveBayes.default}
\alias{naiveBayes.formula}
\alias{print.naiveBayes}
\alias{predict.naiveBayes}
\title{Naive Bayes Classifier}
\description{
  Computes the conditional a-posterior probabilities of a categorical
  class variable given independent predictor variables using
  the Bayes rule.
}
\usage{
\S3method{naiveBayes}{formula}(formula, data, laplace = 0, ..., subset, na.action = na.pass)
\S3method{naiveBayes}{default}(x, y, laplace = 0, ...)


\S3method{predict}{naiveBayes}(object, newdata,
  type = c("class", "raw"), threshold = 0.001, eps = 0, ...)

}
\arguments{
  \item{x}{A numeric matrix, or a data frame of categorical and/or
    numeric variables.}
  \item{y}{Class vector.}
  \item{formula}{A formula of the form \code{class ~ x1 + x2 +
      \dots}. Interactions are not allowed.}
  \item{data}{Either a data frame of predictors (categorical and/or
    numeric) or a contingency table.}
  \item{laplace}{positive double controlling Laplace smoothing. The
    default (0) disables Laplace smoothing.}
  \item{\dots}{Currently not used.}
  \item{subset}{For data given in a data frame, an index vector
    specifying the cases to be used in the
    training sample.  (NOTE: If given, this argument must be
    named.)}
  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is not to count them for the
	  computation of the probability factors. An
          alternative is na.omit, which leads to rejection of cases
          with missing values on any required variable.  (NOTE: If
          given, this argument must be named.)}
  \item{object}{An object of class \code{"naiveBayes"}.}
  \item{newdata}{A dataframe with new predictors (with possibly fewer
    columns than the training data). Note that the column names of
    \code{newdata} are matched against the training data ones.}
  \item{type}{If \code{"raw"}, the conditional a-posterior
    probabilities for each class are returned, and the class with
    maximal probability else.}
  \item{threshold}{Value replacing cells with probabilities within \code{eps} range.}
  \item{eps}{double for specifying an epsilon-range to apply laplace
    smoothing (to replace zero or close-zero probabilities by \code{theshold}.)}
}
\value{
  An object of class \code{"naiveBayes"} including components:

  \item{apriori}{Class distribution for the dependent variable.}
  \item{tables}{A list of tables, one for each predictor variable. For each
    categorical variable a table giving, for each attribute level, the conditional
    probabilities given the target class. For each numeric variable, a
    table giving, for each target class, mean and standard deviation of
    the (sub-)variable.}

}
\details{
  The standard naive Bayes classifier (at least this implementation)
  assumes independence of the predictor
  variables, and Gaussian distribution (given the target class) of
  metric predictors.
  For attributes with missing values, the
  corresponding table entries are omitted for prediction.
  }
\author{David Meyer \email{David.Meyer@R-project.org}. Laplace smoothing
enhancement by Jinghao Xue.}
\examples{
## Categorical data only:
data(HouseVotes84, package = "mlbench")
model <- naiveBayes(Class ~ ., data = HouseVotes84)
predict(model, HouseVotes84[1:10,])
predict(model, HouseVotes84[1:10,], type = "raw")

pred <- predict(model, HouseVotes84)
table(pred, HouseVotes84$Class)

## using laplace smoothing:
model <- naiveBayes(Class ~ ., data = HouseVotes84, laplace = 3)
pred <- predict(model, HouseVotes84[,-1])
table(pred, HouseVotes84$Class)


## Example of using a contingency table:
data(Titanic)
m <- naiveBayes(Survived ~ ., data = Titanic)
m
predict(m, as.data.frame(Titanic))

## Example with metric predictors:
data(iris)
m <- naiveBayes(Species ~ ., data = iris)
## alternatively:
m <- naiveBayes(iris[,-5], iris[,5])
m
table(predict(m, iris), iris[,5])
}
\keyword{classif}
\keyword{category}

