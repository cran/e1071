\name{write.svm}
\alias{write.svm}
\title{Write SVM object to file}
\description{
  This function exports an SVM object (trained by \code{svm}) to a specified file
  in the format that the function 'svm\_load\_model' of libsvm can read.
}
\usage{
write.svm(object, file)
}
\arguments{
  \item{object}{Object of class \code{"svm"}, created by \code{svm}.}
  \item{file}{filename to export the svm object to.}
}
\details{
  This function is useful when SVM models trained in R shall be used in 
  other environments.
}
\author{
  Tomomi TAKASHINA (based on 'predict.svm' by David Meyer)
  \email{t.takashina@computer.org}
}
\seealso{
  \code{\link{svm}}
}
\examples{
data(iris)
attach(iris)

## classification mode
# default with factor response:
model <- svm (Species~., data=iris)

# export SVM object to file
write.svm(model, file="iris-classifier.svm")
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}


