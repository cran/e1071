\name{predict.naiveBayes}
\alias{predict.naiveBayes}
\title{Naive Bayes Classifier}
\description{
  Computes the conditional a-posterior probabilities of a categorical
  class variable given independent predictor variables using
  the Bayes rule.
}
\usage{
predict.naiveBayes(object, newdata, type = c("class", "raw"), threshold = 0.001, ...)
}
\arguments{
  \item{object}{An object of class \code{"naiveBayes"}.}
  \item{newdata}{A dataframe with new predictors.}
  \item{type}{see value.}
  \item{threshold}{Value replacing cells with 0 probabilities.}
  \item{\dots}{Currently not used.}
}
\value{
If \code{type = "raw"}, the conditional a-posterior
probabilities for each class are returned, and the class with
maximal probability else.
}
\details{
  The standard naive Bayes classifier (at least this implementation)
  assumes independence of the predictor
  variables, and gaussian distribution (given the target class) of metric predictors.
  For attributes with missing values, the
  corresponding table entries are omitted for prediction.
}
\author{David Meyer \email{david.meyer@ci.tuwien.ac.at}}
\seealso{\code{\link{naiveBayes}}}
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
\keyword{classif}
\keyword{category}
