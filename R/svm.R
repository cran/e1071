svm <- function (x, ...)
  UseMethod ("svm")

svm.formula <-
function (formula, data = NULL, ..., subset, na.action = na.omit, scale = TRUE)
{
  call <- match.call()
  if (!inherits(formula, "formula")) 
    stop("method is only for formula objects")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$scale <- NULL
  m[[1]] <- as.name("model.frame")
  m$na.action <- na.action
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, response)
  if (length(scale) == 1)
    scale <- rep(scale, ncol(x))
  if (any(scale)) {
    remove <- unique(c(which(labels(Terms) %in% names(attr(x, "contrasts"))),
                       which(!scale)
                       )
                     )
    scale <- !attr(x, "assign") %in% remove
  }
  ret <- svm.default (x, y, scale = scale, ...)
  ret$call <- call
  ret$call[[1]] <- as.name("svm")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) 
    ret$na.action <- attr(m, "na.action")
  class(ret) <- c("svm.formula", class(ret))
  return (ret)
}

svm.default <-
function (x,
          y           = NULL,
          scale       = TRUE,
          type        = NULL,
          kernel      = "radial",
          degree      = 3,
          gamma       = 1 / ncol(as.matrix(x)),
          coef0       = 0,
          cost        = 1,
          nu          = 0.5,
          class.weights = NULL,
          cachesize   = 40,
          tolerance   = 0.001,
          epsilon     = 0.1,
          shrinking   = TRUE,
          cross       = 0,
          probability = FALSE,
          fitted      = TRUE,
          ...,
          subset,
          na.action = na.omit)
{
  sparse  <- inherits(x, "matrix.csr")
  if (sparse) {
    if (!require(SparseM))
      stop("Need SparseM package for handling of sparse structures!")
  }

  ## NULL parameters?
  if(is.null(degree)) stop("`degree' must not be NULL!")
  if(is.null(gamma)) stop("`gamma' must not be NULL!")
  if(is.null(coef0)) stop("`coef0' must not be NULL!")
  if(is.null(cost)) stop("`cost' must not be NULL!")
  if(is.null(nu)) stop("`nu' must not be NULL!")
  if(is.null(epsilon)) stop("`epsilon' must not be NULL!")
  if(is.null(tolerance)) stop("`tolerance' must not be NULL!")
  
  xhold   <- if (fitted) x else NA
  x.scale <- y.scale <- NULL
  formula <- inherits(x, "svm.formula")
  
  ## determine model type
  if (is.null(type)) type <-
    if (is.null(y)) "one-classification"
    else if (is.factor(y)) "C-classification"
    else "eps-regression"

  type <- pmatch(type, c("C-classification",
                         "nu-classification",
                         "one-classification",
                         "eps-regression",
                         "nu-regression"), 99) - 1

  if (type > 10) stop("wrong type specification!")
  
  kernel <- pmatch(kernel, c("linear",
                             "polynomial",
                             "radial",
                             "sigmoid"), 99) - 1

  if (kernel > 10) stop("wrong kernel specification!")

  ## scaling, subsetting, and NA handling
  if (sparse) {
    scale <- rep(FALSE, ncol(x))
    if(!is.null(y)) na.fail(y)
    x <- t(t(x)) ## make shure that col-indices are sorted
  } else {
    x <- as.matrix(x)

    ## subsetting and na-handling for matrices
    if (!formula) {
      if (!missing(subset)) x <- x[subset,]
      if (is.null(y))
        x <- na.action(x)
      else {
        df <- na.action(data.frame(y, x))
        y <- df[,1]
        x <- as.matrix(df[,-1])
      }
    }

    ## scaling
    if (length(scale) == 1)
      scale <- rep(scale, ncol(x))
    if (any(scale)) {
      co <- !apply(x[,scale, drop = FALSE], 2, var)
      if (any(co)) {
        scale <- rep(FALSE, ncol(x))
        warning(paste("Variable(s)",
                      paste("`",colnames(x[,scale, drop = FALSE])[co],
                            "'", sep="", collapse=" and "),
                      "constant. Cannot scale data.")
                )
      } else {
        xtmp <- scale(x[,scale])
        x[,scale] <- xtmp
        x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
        if (is.numeric(y) && (type > 2)) {
          y <- scale(y)
          y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
          y <- as.vector(y)
        }
      }
    }
  }

  ## further parameter checks
  nr <- nrow(x)
  if (cross > nr)
    stop("`cross' cannot exceed the number of observations!")

  if (!is.vector(y) && !is.factor (y) && type != 2) stop("y must be a vector or a factor.")
  if (type != 2 && length(y) != nr) stop("x and y don't match.")

  if (cachesize < 0.1) cachesize <- 0.1

  if (type > 2 && !is.numeric(y))
    stop("Need numeric dependent variable for regression.")

  lev <- NULL
  weightlabels <- NULL
  # in case of classification: transform factors into integers
  if (type == 2) # one class classification --> set dummy
    y <- 1
  else
    if (is.factor(y)) {
      lev <- levels(y)
      y <- as.integer(y)
      if (!is.null(class.weights)) {
        if (is.null(names(class.weights)))
          stop("Weights have to be specified along with their according level names !")
        weightlabels <- match (names(class.weights), lev)
        if (any(is.na(weightlabels)))
          stop("At least one level name is missing or misspelled.")
      }
    } else {
      if (type < 3 && any(as.integer(y) != y))
        stop("dependent variable has to be of factor or integer type for classification mode.")
      lev <- unique(y)
   }

  nclass <- 2
  if (type < 2) nclass <- length(lev)

  if (type > 1 && length(class.weights) > 0) {
    class.weights <- NULL
    warning("`class.weights' are set to NULL for regression mode. For classification, use a _factor_ for `y', or specify the correct `type' argument.")
  }

  cret <- .C ("svmtrain",
              # data
              as.double  (if (sparse) x@ra else t(x)),
              as.integer (nr), as.integer(ncol(x)),
              as.double  (y),

              # sparse index info
              as.integer (if (sparse) x@ia else 0),
              as.integer (if (sparse) x@ja else 0), 
              
              # parameters
              as.integer (type),
              as.integer (kernel),
              as.double  (degree),
              as.double  (gamma),
              as.double  (coef0),
              as.double  (cost),
              as.double  (nu),
              as.integer (weightlabels),
              as.double  (class.weights),
              as.integer (length (class.weights)),
              as.double  (cachesize),
              as.double  (tolerance),
              as.double  (epsilon),
              as.integer (shrinking),
              as.integer (cross),
              as.integer (sparse),
              as.integer (probability),

              # results
              nclasses = integer  (1), 
              nr       = integer  (1), # nr of support vectors
              index    = integer  (nr),
              labels   = integer  (nclass),
              nSV      = integer  (nr),
              rho      = double   (nclass * (nclass - 1) / 2),
              coefs    = double   (nr * (nclass - 1)),
              sigma    = double   (1),
              probA    = double   (nclass * (nclass - 1) / 2),
              probB    = double   (nclass * (nclass - 1) / 2),
              
              cresults = double   (cross),
              ctotal1  = double   (1),
              ctotal2  = double   (1),
              error    = character(1),

              PACKAGE = "e1071")

  if (nchar(cret$error))
    stop(paste(cret$error, "!", sep=""))

  ret <- list (
               call     = match.call(),
               type     = type,
               kernel   = kernel,
               cost     = cost,
               degree   = degree,
               gamma    = gamma,
               coef0    = coef0,
               nu       = nu,
               epsilon  = epsilon,
               sparse   = sparse,
               scaled   = scale,
               x.scale  = x.scale,
               y.scale  = y.scale,
               
               nclasses = cret$nclasses,            #number of classes
               levels   = lev,
               tot.nSV  = cret$nr,                  #total number of sv
               nSV      = cret$nSV[1:cret$nclasses],#number of SV in diff. classes
               labels   = cret$label[1:cret$nclasses],#labels of the SVs.
               SV       = t(t(x[cret$index[1:cret$nr],])), #copy of SV
               index    = cret$index[1:cret$nr],     #indexes of sv in x
               #constants in decision functions
               rho      = cret$rho[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
               #probabilites
               compprob = probability,
               probA    = if (!probability) NULL else
                             cret$probA[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
               probB    = if (!probability) NULL else
                             cret$probB[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
               sigma    = if (probability) cret$sigma else NULL,
               #coefficiants of sv
               coefs    = if (cret$nr == 0) NULL else
                              t(matrix(cret$coefs[1:((cret$nclasses - 1) * cret$nr)],
                                       nrow = cret$nclasses - 1,
                                       byrow = TRUE))
              )

  # cross-validation-results
  if (cross > 0)    
    if (type > 2) {
      scale.factor     <- if (any(scale)) crossprod(y.scale$"scaled:scale") else 1;
      ret$MSE          <- cret$cresults * scale.factor;
      ret$tot.MSE      <- cret$ctotal1  * scale.factor;
      ret$scorrcoeff   <- cret$ctotal2;
    } else {
      ret$accuracies   <- cret$cresults;
      ret$tot.accuracy <- cret$ctotal1;
    }

  class (ret) <- "svm"
  ret$fitted  <- if (fitted) predict(ret, xhold) else NA
  ret
} 

predict.svm <- function (object, newdata,
                         decision.values = FALSE,
                         probability = FALSE,
                         ...,
                         na.action = na.omit) {
  if (missing(newdata))
    return(fitted(object))

  sparse <- inherits(newdata, "matrix.csr")
  if (object$sparse || sparse) {
    if (!require(SparseM))
      stop("Need SparseM package for handling of sparse structures!")
  }

  if (is.vector(newdata) || sparse) newdata <- t(t(newdata))
  if (!object$sparse) {
    if (inherits(object, "svm.formula")) {
      if(is.null(colnames(newdata)))
        colnames(newdata) <- colnames(object$SV)
      newdata <- model.matrix(delete.response(terms(object)),
                              as.data.frame(newdata), na.action = na.action)
    } else if (!sparse) newdata <- na.action(as.matrix(newdata))
  }

  if (any(object$scaled))
    newdata[,object$scaled] <-
      scale(newdata[,object$scaled, drop = FALSE],
            center = object$x.scale$"scaled:center",
            scale  = object$x.scale$"scaled:scale"
            )

  if (ncol(object$SV) != ncol(newdata)) stop ("test data does not match model !")

  ret <- .C ("svmpredict",
             as.integer (decision.values),
             as.integer (probability),
             
             # model
             as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
             as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
             as.integer (if (object$sparse) object$SV@ia else 0),
             as.integer (if (object$sparse) object$SV@ja else 0),
             as.double  (as.vector(object$coefs)),
             as.double  (object$rho),
             as.integer (object$compprob),
             as.double  (object$probA),
             as.double  (object$probB),
             as.integer (object$nclasses),
             as.integer (object$tot.nSV),
             as.integer (object$labels),
             as.integer (object$nSV),
             as.integer (object$sparse),
             
             # parameter
             as.integer (object$type),
             as.integer (object$kernel),
             as.double  (object$degree),
             as.double  (object$gamma),
             as.double  (object$coef0),

             # test matrix
             as.double  (if (sparse) newdata@ra else t(newdata)),
             as.integer (nrow(newdata)),
             as.integer (if (sparse) newdata@ia else 0),
             as.integer (if (sparse) newdata@ja else 0),
             as.integer (sparse),
             
             # decision-values
             ret = double(nrow(newdata)),
             dec = double(nrow(newdata) * object$nclasses * (object$nclasses - 1) / 2),
             prob = double(nrow(newdata) * object$nclasses),

             PACKAGE = "e1071"
            )
  
  ret2 <- if (is.character(object$levels))
    # classification: return factors
    factor (object$levels[ret$ret], levels = object$levels)
  else if (object$type == 2)
    # one-class-classification: return TRUE/FALSE
    ret$ret == 1 
  else if (any(object$scaled))
    # return raw values, possibly scaled back
    ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
  else
    ret$ret

  if (decision.values) {
    colns = c()
    for (i in 1:(object$nclasses - 1))
      for (j in (i + 1):object$nclasses)
        colns <- c(colns, paste(object$levels[i],"/",object$levels[j], sep = ""))
    attr(ret2, "decision.values") <- matrix(ret$dec, nrow = nrow(newdata), byrow = TRUE)
    colnames(attr(ret2, "decision.values")) <- colns
  }

  if (probability && object$type < 2) {
    attr(ret2, "probabilities") <- matrix(ret$prob, nrow = nrow(newdata), byrow = TRUE)
    colnames(attr(ret2, "probabilities")) <- object$levels
  }

  ret2
}

print.svm <- function (x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("Parameters:\n")
  cat("   SVM-Type: ", c("C-classification",
                         "nu-classification",
                         "one-classification",
                         "eps-regression",
                         "nu-regression")[x$type+1], "\n")
  cat(" SVM-Kernel: ", c("linear",
                         "polynomial",
                         "radial",
                         "sigmoid")[x$kernel+1], "\n")
  if (x$type==0 || x$type==3 || x$type==4)
    cat("       cost: ", x$cost, "\n")
  if (x$kernel==1)
    cat("     degree: ", x$degree, "\n")
  cat("      gamma: ", x$gamma, "\n")
  if (x$kernel==1 || x$kernel==3)
    cat("     coef.0: ", x$coef0, "\n")
  if (x$type==1 || x$type==2 || x$type==4)
    cat("         nu: ", x$nu, "\n")
  if (x$type==3) {
    cat("    epsilon: ", x$epsilon, "\n\n")
    if (x$compprob)
      cat("Sigma: ", x$sigma, "\n\n")
  }
  
  cat("\nNumber of Support Vectors: ", x$tot.nSV)
  cat("\n\n")
  
}

summary.svm <- function(object, ...)
  structure(object, class="summary.svm")

print.summary.svm <- function (x, ...) {
  print.svm(x)
  if (x$type<2) {
    cat(" (", x$nSV, ")\n\n")
    cat("\nNumber of Classes: ", x$nclasses, "\n\n")
    cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels)
  }
  cat("\n\n")
  if (x$type==2) cat("\nNumber of Classes: 1\n\n\n")

  if ("MSE" %in% names(x)) {
    cat(length (x$MSE), "-fold cross-validation on training data:\n\n", sep="")
    cat("Total Mean Squared Error:", x$tot.MSE, "\n")
    cat("Squared Correlation Coefficient:", x$scorrcoef, "\n")
    cat("Mean Squared Errors:\n", x$MSE, "\n\n")
  }
  if ("accuracies" %in% names(x)) {
    cat(length (x$accuracies), "-fold cross-validation on training data:\n\n", sep="")
    cat("Total Accuracy:", x$tot.accuracy, "\n")
    cat("Single Accuracies:\n", x$accuracies, "\n\n")
  }
  cat("\n\n")
}

scale.data.frame <- function(x, center = TRUE, scale = TRUE) {
  i <- sapply(x, is.numeric)
  if (ncol(x[, i, drop = FALSE])) {
    x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
    if(center || !is.logical(center))
      attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
    if(scale || !is.logical(scale))
      attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
  }
  x
}

plot.svm <- function(x, data, formula = NULL, fill = TRUE,
                     grid = 50, slice = list(), ...) {
  if (x$type < 3) {
    if (is.null(formula) && ncol(data) == 3) {
      formula <- formula(delete.response(terms(x)))
      formula[2:3] <- formula[[2]][2:3]
    }
    if(is.null(formula)) stop("missing formula.")
    if (fill) {
      sub <- model.frame(formula, data)
      xr <- seq(min(sub[,2]), max(sub[,2]), length = grid)
      yr <- seq(min(sub[,1]), max(sub[,1]), length = grid)
      l <- length(slice)
      if (l < ncol(data) - 3) {
        slnames <- names(slice)
        slice <- c(slice, rep(list(0), ncol(data) - 3 - l))
        names <- labels(delete.response(terms(x)))
        names(slice) <- c(slnames, names[! names %in% c(colnames(sub), slnames)])
      }
      lis <- c(list(yr), list(xr), slice)
      names(lis)[1:2] <- colnames(sub)
      new <- expand.grid(lis)[,labels(terms(x))]
      preds <- predict(x, new)
      filled.contour(xr, yr, matrix(as.numeric(preds), nr = length(xr), byrow=TRUE),
                     plot.axes = {
                       axis(1)
                       axis(2)
                       colors <- as.numeric(model.response(model.frame(x, data)))
                       points(formula, data = data[-x$index,], col = colors[-x$index])
                       points(formula, data = data[x$index,], pch = "x", col = colors[x$index])
                     },
                     levels = 1:(length(unique(as.numeric(preds)))+1),
                     key.axes = axis(4,
                       1:length(unique(as.numeric(preds)))+0.5,
                       labels = levels(preds)[unique(preds)], las = 3
                       ),
                     plot.title = title(main = "SVM classification plot",
                       xlab = names(lis)[2], ylab = names(lis)[1]),
                     ...
                     )
    } else {
      plot(formula, data = data, type = "n", ...)
      colors <- as.numeric(model.response(model.frame(m, data)))
      points(formula, data = data[-x$index,], col = colors[-x$index])
      points(formula, data = data[x$index,], pch = "x", col = colors[x$index])
    }
  }
}

write.svm <- function (object, svm.file="Rdata.svm", scale.file = "Rdata.scale") {

  ret <- .C ("svmwrite",
             # model
             as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
             as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
             as.integer (if (object$sparse) object$SV@ia else 0),
             as.integer (if (object$sparse) object$SV@ja else 0),
             as.double  (as.vector(object$coefs)),
             as.double  (object$rho),
             as.double  (object$probA),
             as.double  (object$probB),
             as.integer (object$nclasses),
             as.integer (object$tot.nSV),
             as.integer (object$labels),
             as.integer (object$nSV),
             as.integer (object$sparse),
             
             # parameter
             as.integer (object$type),
             as.integer (object$kernel),
             as.double  (object$degree),
             as.double  (object$gamma),
             as.double  (object$coef0),

             # filename
             as.character(svm.file),

             PACKAGE = "e1071"
            )$ret

   write.table(data.frame(center = object$x.scale$"scaled:center",
                          scale  = object$x.scale$"scaled:scale"),
               file=scale.file, col.names=FALSE, row.names=FALSE)
}
