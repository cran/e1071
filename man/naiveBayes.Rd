\name{naiveBayes}
\alias{naiveBayes}
\alias{naiveBayes.default}
\alias{naiveBayes.formula}
\alias{print.naiveBayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Naive Bayes Classifier}
\description{
  Computes the conditional a-posterior probabilities of a categorical
  class variable given independent predictor variables using
  the Bayes rule.
}
\usage{
\method{naiveBayes}{formula}(formula, data, ..., subset, na.action = na.pass)
\method{naiveBayes}{default}(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric matrix, or a data frame of categorical and/or
    numeric variables.}
  \item{y}{Class vector.}
  \item{formula}{A formula of the form \code{class ~ x1 + x2 +
      \dots}. Interactions are not allowed.}
  \item{data}{Either a data frame of predictors (caegorical and/or
    numeric) or a contingency table.}
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
  variables, and gaussian distribution (given the target class) of metric predictors.
  }
\author{David Meyer \email{david.meyer@ci.tuwien.ac.at}}
\seealso{\code{\link{predict.naiveBayes}}}
\examples{
## Categorical data only:
data(HouseVotes84)
model <- naiveBayes(Class ~ ., data = HouseVotes84)
predict(model, HouseVotes84[1:10,-1])
predict(model, HouseVotes84[1:10,-1], type = "raw")

pred <- predict(model, HouseVotes84[,-1])
table(pred, HouseVotes84$Class)

## Example of using a contingency table:
data(Titanic)
m <- naiveBayes(Survived ~ ., data = Titanic)
m
predict(m, as.data.frame(Titanic)[,1:3])

## Example with metric predictors:
data(iris)
m <- naiveBayes(Species ~ ., data = iris)
## alternatively:
m <- naiveBayes(iris[,-5], iris[,5])
m
table(predict(m, iris[,-5]), iris[,5])
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{category}% at least one, from doc/KEYWORDS

