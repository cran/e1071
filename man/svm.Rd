\name{svm}
\alias{svm}
\alias{svm.default}
\alias{svm.formula}
\alias{summary.svm}
\alias{print.summary.svm}
\alias{print.svm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Support Vector Machines}
\description{
\code{svm} is used to train a support vector machine. It can be used to carry
out general regression and classification (of nu and epsilon-type), as
well as density-estimation. A formula interface is provided.
}
\usage{
\method{svm}{formula}(formula, data = list(), subset, na.action = na.fail,
...)
\method{svm}{default}(x, y = NULL, scale = TRUE, subset, na.action =
na.fail, type = NULL, kernel = "radial", degree = 3, gamma = 1 / dim(x)[2],
coef0 = 0, cost = 1, nu = 0.5, class.weights = NULL, cachesize = 40,
tolerance = 0.001, epsilon = 0.5, shrinking = TRUE, cross = 0, fitted = TRUE, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. Note, that an
          intercept is always included, whether given in the formula or
          not.}
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `svm' is called from.}
  \item{x}{a data matrix, a vector, or a sparse matrix (object of class
    \code{\link[SparseM]{matrix.csr}} as provided by the package \code{SparseM}).}
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    regression).}
  \item{scale}{Per default, data are scaled internally (both \code{x} and \code{y}
    variables) to zero mean and unit variance. The center and scale
    values are returned and used for later predictions.}
  \item{type}{\code{svm} can be used as a classification
    machine, as a regresson machine or a density estimator. Depending of whether \code{y} is
    a factor or not, the default setting for \code{type} is \code{C-classification} or \code{eps-regression}, respectively, but may be overwritten by setting an explicit value.\cr
    Valid options are:
    \itemize{
      \item \code{C-classification}
      \item \code{nu-classification}
      \item \code{one-classification} (for density estimation)
      \item \code{eps-regression}
      \item \code{nu-regression}
    }
  }
  \item{kernel}{the kernel used in training and predicting. You
    might consider changing some of the following parameters, depending
    on the kernel type.\cr
    \describe{
      \item{linear:}{\eqn{u'v}{u'*v}}
      \item{polynomial:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
      \item{radial basis:}{\eqn{e^(-\gamma |u-v|^2)}{exp(-gamma*|u-v|^2)}}
      \item{sigmoid:}{\eqn{tanh(\gamma u'v + coef0)}{tanh(gamma*u'*v + coef0)}}
      }
    }
  \item{degree}{parameter needed for kernel of type \code{polynomial} (default: 3)}
  \item{gamma}{parameter needed for all kernels except \code{linear}
    (default: 1/(data dimension))}
  \item{coef0}{parameter needed for kernels of type \code{polynomial}
    and \code{sigmoid} (default: 0)}
  \item{cost}{cost of constraints violation (default: 1)---it is the
    `C'-constant of the regularization term in the Lagrange formulation.}
  \item{nu}{parameter needed for \code{nu-classification} and \code{one-classification}}
  \item{class.weights}{a named vector of weights for the different
    classes, used for asymetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}
  \item{cachesize}{cache memory in MB (default 40)}
  \item{tolerance}{tolerance of termination criterion (default: 0.001)}
  \item{epsilon}{epsilon in the insensitive-loss function (default: 0.5)}
  \item{shrinking}{option whether to use the shrinking-heuristics
    (default: \code{TRUE})}
  \item{cross}{if a integer value k>0 is specified, a k-fold cross
    validation on the training data is performed to assess the quality
    of the model: the accuracy rate for classification and the Mean
    Sqared Error for regression}
  \item{fitted}{indicates whether the fitted values should be computed
    and included in the model or not (default: \code{TRUE})}
  \item{\dots}{additional parameters for the low level fitting function
    \code{svm.default}}
  \item{subset}{An index vector specifying the cases to be used in the
          training sample.  (NOTE: If given, this argument must be
          named.)}
  \item{na.action}{A function to specify the action to be taken if `NA's are
          found. The default action is for the procedure to fail.  An
          alternative is `na.omit', which leads to rejection of cases
          with missing values on any required variable.  (NOTE: If
          given, this argument must be named.)}	
}

\value{
  An object of class \code{"svm"} containing the fitted model, including:
  \item{SV}{the resulting support vectors (possibly scaled)}
  \item{index}{the index of the resulting support vectors in the data matrix}
  \item{coefs}{the corresponding coefficients times the training labels}
  \item{rho}{the negative intercept}
  (Use \code{summary} and \code{print} to get some output).
}
\details{For multiclass-classification with k levels, k>2, \code{libsvm} uses the
  `one-against-one'-approach, in which k(k-1)/2 binary classifiers are
  trained; the appropriate class is found by a voting scheme.\cr
\code{libsvm} internally uses a sparse data representation, which is 
also high-level supported by the package \code{SparseM}.\cr
If the predictor variables include factors, the formula interface must be used to get a
correct model matrix. \code{plot.svm} allows a simple graphical
visualization of classification models.}
\note{Data are scaled internally, usually yielding better results.}
\references{
  \itemize{
    \item
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}

     \item
      Exact formulations of models, algorithms, etc. can be found in the
      document:\cr
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.ps.gz}
    
    \item 
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{Libsvm: Introduction and Benchmarks}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/q2.ps.gz}
    
  }
}
\author{
  David Meyer (based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin)\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{
  \code{\link{predict.svm}}
  \code{\link{plot.svm}}
  \code{\link[SparseM]{matrix.csr}} (in package `SparseM')
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
newdata
predict (m, newdata)

# visualization:
plot (as.matrix(X))
points (as.matrix(X)[m$index,], col=2)
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}









