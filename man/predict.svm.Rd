\name{predict.svm}
\alias{predict.svm}
\title{Predict method for Support Vector Machines}
\description{
  This function predicts values based upon a model trained by \code{svm}.
}
\usage{
\method{predict}{svm}(object, newdata, ..., na.action = na.omit)
}

\arguments{
  \item{object}{Object of class \code{"svm"}, created by \code{svm}.}
  \item{newdata}{A matrix containing the new input data. A vector will
    be transformed to a n x 1 matrix.}
  \item{na.action}{A function to specify the action to be taken if `NA's are
    found. The default action is \code{na.omit}, which leads to rejection of cases
    with missing values on any required variable. An alternative
    is \code{na.fail}, which causes an error if \code{NA} cases
    are found. (NOTE: If given, this argument must be named.)}	
  \item{\dots}{Currently not used.}
}
\value{
  The predicted value (for classification: the label, for density
  estimation: \code{TRUE} or \code{FALSE}).
}
\note{
  If the training set was scaled by \code{svm} (done by default), the
  new data is scaled accordingly using scale and center of
  the training data.
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
attach(iris)

## classification mode
# default with factor response:
model <- svm (Species~., data=iris)

# alternatively the traditional interface:
x <- subset (iris, select = -Species)
y <- Species
model <- svm (x, y) 

print (model)
summary (model)

# test with train data
pred <- predict (model, x)
# (same as:)
pred <- predict (model)

# Check accuracy:
table (pred,y)

## try regression mode on two dimensions

# create data
x <- seq (0.1,5,by=0.05)
y <- log(x) + rnorm (x, sd=0.2)

# estimate model and predict input values
m   <- svm (x,y)
new <- predict (m,x)

# visualize
plot   (x,y)
points (x, log(x), col=2)
points (x, new, col=4)

## density-estimation

# create 2-dim. normal with rho=0:
X <- data.frame (a=rnorm (1000), b=rnorm (1000))
attach (X)

# traditional way:
m <- svm (X, gamma=0.1)

# formula interface:
m <- svm (~., data=X, gamma=0.1)
# or:
m <- svm (~a+b, gamma=0.1)

# test:
newdata <- data.frame(a=c(0,4), b=c(0,4))
predict (m, newdata)
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}

