\name{write.svm} 
\alias{write.svm}
\title{Write SVM object to file}
\description{
  This function exports an SVM object (trained by \code{svm}) to two
  specified files. One is in the format that the
  function 'svm\_load\_model' of libsvm can read. The other is for
  scaling data, containing a data with centers and scales for all variables.
}
\usage{
write.svm(object, svm.file, scale.file)
}
\arguments{
  \item{object}{Object of class \code{"svm"}, created by \code{svm}.}
  \item{svm.file}{filename to export the svm object to.}
  \item{scale.file}{filename to export the scaling data of the svm object to.}
}
\details{
  This function is useful when SVM models trained in R shall be used in 
  other environments. The SVM model is saved in the standard format of
  libsvm. The scaling data is written to a separate file because scaling
  data is not included in the standard format of libsvm. The format
  of the scaling data file is a n times 2 matrix: the n-th row
  corresponds to the n-th dimension of the data, the colums being formed
  of the corresponding mean and scale.
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
write.svm(model, svm.file = "iris-classifier.svm", scale.file = "iris-classifier.scale")

# read scale file
# the n-th row is corresponding to n-th dimension. The 1st column contains the
# center value, the 2nd column is the scale value.

read.table("iris-classifier.scale")

}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}

