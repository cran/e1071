\name{svm}
\alias{svm}
\alias{svm.default}
\alias{svm.formula}
\alias{summary.svm}
\alias{print.svm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Support Vector Machines}
\description{
\code{svm} is used to train a support vector machine. It can be used to carry
out general regression and classification (of nu and epsilon-type), as
well as density-estimation. A formula interface is provided.
}
\usage{
\method{svm}{formula}(formula, data = NULL, ...)
\method{svm}{default}(x, y=NULL, type=NULL, kernel="radial", degree=3, gamma=1/dim(x)[2],
coef0=0, cost=1, nu=0.5, class.weights=NULL, cachesize=40, tolerance=0.001, epsilon=0.5,
shrinking=TRUE, cross=0, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. Note, that an
          intercept is always included, whether given in the formula or
          not.}
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `svm' is called from.}
  \item{x}{a data matrix or a vector.}
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
  a factor (for classification tasks) or a numeric vector (for regression).}
  \item{type}{\code{svm} can be used as a classification
    machine, as a regresson machine or a density estimator. Depending of whether \code{y} is
    a factor or not, the default setting for \code{svm.type} is \code{C-classification} or \code{eps-regression}, respectively, but may be overwritten by setting an explicit value.\cr
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
  \item{cost}{cost of constraints violation. (default: 1)}
  \item{nu}{parameter needed for \code{nu-classification} and \code{one-classification}}
  \item{class.weights}{a named vector of weights for the different
    classes, used for asymetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}
  \item{cachesize}{cache memory in MB. (default 40)}
  \item{tolerance}{tolerance of termination criterion (default: 0.001)}
  \item{epsilon}{epsilon in the insensitive-loss function (default: 0.5)}
  \item{shrinking}{option whether to use the shrinking-heuristics
    (default: TRUE)}
  \item{cross}{if a integer value k>0 is specified, a k-fold cross
    validation on the training data is performed to assess the quality
    of the model: the accuracy rate for classification and the Mean
    Sqared Error for regression}
  \item{\dots}{additional parameters for the low level fitting function
    \code{svm.default}.}
}

\value{
  An object of class \code{"svm"} containing the fitted model, especially:
  \item{sv}{the resulting support vectors}
  \item{index}{the index of the resulting support vectors in the data matrix}
  \item{coefs}{the corresponding coefficiants}
  (Use \code{summary} and \code{print} to get some output).
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
  David Meyer (based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin)\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{
  \code{\link{predict.svm}}
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
m <- svm (X)

# formula interface:
m <- svm (~a+b)
# or:
m <- svm (~., data=X)

# test:
predict (m, t(c(0,0)))
predict (m, t(c(4,4)))

# visualization:
plot (X)
points (X[m$index,], col=2)
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}









