\name{predict.svm}
\alias{predict.svm}
\title{Predict method for Support Vector Machines}
\description{
  This function predicts values based upon a model trained by \code{svm}.
}
\usage{
predict.svm(model, x, type="class")
}

\arguments{
  \item{model}{object of class \code{svm}, created by \code{svm}.}
  \item{x}{A matrix containing the input data.}
  \item{type}{If \code{raw} is supplied, the estimated values are returned. 
    With \code{class}, \code{predict} returns the first level of the
    training response for negative (non-zero) values, and  the second
    level otherwise. If \code{class} is used on numerical data, a factor
    with labels "-1" and  "1" is created.}
  
}
\value{
  The predicted value.
}
\references{
  \itemize{
    \item
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM 2.0: Solving Different Support Vector Formulations.}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm2.ps.gz}
    
    \item 
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{Libsvm: Introduction and Benchmarks}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/q2.ps.gz}
    
  }
}
\author{
  David Meyer (based on C++-code by Chih-Chung Chang and Chih-Jen Lin)\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{
  \code{\link{svm}}
}
\examples{
data(iris)
# amputate data to two factors
iris.sub <- subset(iris, Species != "virginica")

# get independent vars
x <- subset (iris.sub, select = -Species)

# get responses
y <- iris.sub[,"Species"]

# coercion needed for correct factor levels
y <- as.factor(as.character(y))

# default with factor response: classification mode
model <- svm (x, y)
print (model)
summary (model)

# test with train data
predict (model, x)

# try regression mode on two dimensions in linear mode
model <- svm (x[,"Petal.Length"], x[,"Petal.Width"],
svm.type="regression", kernel.type="linear")
print (model)


pred <- predict (model,x[,"Petal.Length"])

par (mfcol=c(1,2))
plot(x[,"Petal.Length"],x[,"Petal.Width"])
plot(x[,"Petal.Length"],pred)

}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}
