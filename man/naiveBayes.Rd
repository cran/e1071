\name{naiveBayes}
\alias{naiveBayes}
\alias{print.naiveBayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Naive Bayes Classifier}
\description{
  Computes the conditional a-posterior probabilities of a categorical
  class variable given independent categorical predictor variables using
  the Bayes rule.
}
\usage{
naiveBayes(formula, data, ..., subset, na.action = na.pass)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula of the form \code{class ~ x1 + x2 +
      \dots}. Interactions are not allowed.}
  \item{data}{Either a data frame of factors or a contingency table.}
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
  \item{tables}{A list of probability tables,
    one for each predictor variable, giving, for each attribute level, the conditional
    probabilities given the predictor classes.}

}
\author{David Meyer \email{david.meyer@ci.tuwien.ac.at}}
\seealso{\code{\link{predict.naiveBayes}}}
\examples{
data(HouseVotes84)
model <- naiveBayes(Class ~ ., data = HouseVotes84)
predict(model, HouseVotes84[1:10,-1])
predict(model, HouseVotes84[1:10,-1], type = "raw")

pred <- predict(model, HouseVotes84[,-1])
table(pred, HouseVotes84$Class)

data(Titanic)
m <- naiveBayes(Survived ~ ., data = Titanic)
m
predict(m, as.data.frame(Titanic)[,1:3])

}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{category}% at least one, from doc/KEYWORDS

