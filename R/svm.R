svm <- function (x, ...)
  UseMethod ("svm")

svm.formula <-
function (formula, data=NULL, ...)
{
  call <- match.call()
  if (!inherits(formula, "formula")) 
    stop("method is only for formula objects")
  m <- match.call(expand = FALSE)
  if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, response)
  ret <- svm.default (x, y, ...)
  attr (ret, "call") <- call
  return (ret)
}

svm.default <-
function (x,
          y         = NULL,
          type      = NULL,
          kernel    = "radial",
          degree    = 3,
          gamma     = 1/ncol(as.matrix(x)),
          coef0     = 0,
          cost      = 1,
          nu        = 0.5,
          class.weights = NULL,
          cachesize = 40,
          tolerance = 0.001,
          epsilon   = 0.5,
          shrinking = TRUE,
          cross     = 0,
          ...)
{
  if (is.vector(x))
    x <- t(t(x))
  else
    x <- as.matrix(x)

  if (is.null (type)) type <-
    if (is.null(y)) "one-classification"
    else if (is.factor(y)) "C-classification"
    else "eps-regression"

  type <- pmatch (type, c("C-classification",
                          "nu-classification",
                          "one-classification",
                          "eps-regression",
                          "nu-regression"),1) - 1
  
  kernel <- pmatch (kernel, c("linear",
                              "polynomial",
                              "radial",
                              "sigmoid"),3) - 1

  if (!is.vector(y) && !is.factor (y) && !(type==2)) stop ("y must be a vector or a factor.")
  if ((type !=2) && length(y) != nrow(x)) stop ("x and y don't match.")

  if (cachesize < 0.1) cachesize <- 0.1
  
  lev <- NULL
  weightlabels <- NULL
  # in case of classification: transform factors into integers
  if (type == 2) # one class classification --> set dummy
    y <- 1
  else
    if (is.factor(y)) {
      lev <- levels (y)
      y <- codes (y)
      if (!is.null(class.weights)) {
        if (is.null(names (class.weights)))
          stop ("Weights have to be specified along with their according level names !")
        weightlabels <- match (names(class.weights),lev)
        if (any(is.na(weightlabels)))
          stop ("At least one level name is missing or misspelled.")
      }
    } else if (type < 3 && !is.integer (y))
      stop ("dependent variable has to be of factor or integer type for classification mode.")

  nclass <- 2
  if (type < 2) nclass <- length (lev)
  
  cret <- .C ("svmtrain",
              # parameters
              as.double  (t(x)), as.integer(nrow (x)), as.integer(ncol(x)),
              as.double  (y),
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

              # results
              nclasses = integer (1), 
              nr       = integer (1), # nr of support vectors
              index    = integer (nrow(x)),
              labels   = integer (nclass),
              nSV      = integer (nrow(x)),
              rho      = double  (nclass*(nclass-1)/2),
              coefs    = double  (nrow(x)*(nclass-1)),
              
              cresults = double  (cross),
              ctotal1  = double  (1),
              ctotal2  = double  (1)
             )

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
               
               nclasses = cret$nclasses,            #number of classes
               levels   = lev,
               tot.nSV  = cret$nr,                  #total number of sv
               nSV      = cret$nSV[1:cret$nclasses],#number of SV in diff. classes
               labels   = cret$label[1:cret$nclasses],#labels of the SVs.
               SV       = t(t(x[cret$index,])),  #copy of sv
               index    = cret$index,            #indexes of sv in x
               #constants in decision functions
               rho      = cret$rho[1:(cret$nclasses*(cret$nclasses-1)/2)],
               #coefficiants of sv
               coefs    = if (cret$nr==0) NULL else
                              t(matrix(cret$coefs[1:((cret$nclasses-1)*cret$nr)],
                                       nrow=cret$nclasses-1,
                                       byrow=TRUE))
              )

  # cross-validation-results
  if (cross > 0)    
    if (type > 2) {
      ret$MSE <- cret$cresults;
      ret$tot.MSE <- cret$ctotal1;
      ret$scorrcoeff <- cret$ctotal2;
    } else {
      ret$accuracies <- cret$cresults;
      ret$tot.accuracy <- cret$ctotal1;
    }

  class (ret) <- "svm"
  ret
} 

predict.svm <- function (object, newdata, ...) {
  if (is.vector (newdata))
    newdata <- t(t(newdata))
  else
    newdata <- as.matrix(newdata)
    
  if (ncol(newdata) != ncol(object$SV)) stop ("test vector does not match model !")

  ret <- .C ("svmpredict",
             #model
             as.double  (t(object$SV)),
             as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
             as.double  (as.vector (object$coefs)),
             as.double  (object$rho),
             as.integer (object$nclasses),
             as.integer (object$tot.nSV),
             as.integer (object$labels),
             as.integer (object$nSV),
             
             #parameter
             as.integer (object$type),
             as.integer (object$kernel),
             as.double  (object$degree),
             as.double  (object$gamma),
             as.double  (object$coef0),

             #test matrix
             as.double (t(newdata)), as.integer (nrow(newdata)),
             
             #decision-values
             ret = double  (nrow(newdata))
            )$ret

  if (!is.null(object$levels))
    #classification: return factors
    as.factor (object$levels[ret])
  else
    if (object$type == 2)
      #one-class-classification: return TRUE/FALSE
      ret == 1 
    else
      #else: return raw values
      ret
}

print.svm <- function (x, ...) {
  cat ("\nCall:\n",deparse (x$call),"\n\n")
  cat ("Parameters:\n")
  cat ("   SVM-Type: ",c("C-classification",
                         "nu-classification",
                         "one-classification",
                         "eps-regression",
                         "nu-regression")[x$type+1],"\n")
  cat (" SVM-Kernel: ",c("linear",
                         "polynomial",
                         "radial",
                         "sigmoid")[x$kernel+1],"\n")
  cat ("       cost: ",x$cost,"\n")
  cat ("     degree: ",x$degree,"\n")
  cat ("      gamma: ",x$gamma,"\n")
  cat ("     coef.0: ",x$coef0,"\n")
  cat ("         nu: ",x$nu,"\n")
  cat ("    epsilon: ",x$epsilon,"\n")
  cat ("       cost: ",x$cost,"\n\n")
  
  cat ("\nNumber of Support Vectors: ",x$tot.nSV," (",x$nSV,")\n\n")
  if (x$type<2) {
    cat ("\nNumber of Classes: ",x$nclasses,"\n\n")
    cat ("Levels:\n",x$lev,"\n\n")
  }
  if (x$type==2) cat ("\nNumber of Classes: 1\n\n\n")
  cat ("Rho:\n",x$rho,"\n\n")

  if ("MSE" %in% names(x)) {
    cat (length (x$MSE),"-fold cross-validation on training data:\n\n",sep="")
    cat ("Total Mean Squared Error:",x$tot.MSE,"\n")
    cat ("Squared Correlation Coefficient:",x$scorrcoef,"\n")
    cat ("Mean Squared Errors:\n",x$MSE,"\n\n")
  }
  if ("accuracies" %in% names(x)) {
    cat (length (x$accuracies),"-fold cross-validation on training data:\n\n",sep="")
    cat ("Total Accuracy:",x$tot.accuracy,"\n")
    cat ("Single Accuracies:\n",x$accuracies,"\n\n")
  }
}

summary.svm <- function (object, ...) {
  print (object)
  cat ("Support Vectors:\n")
  print (object$SV)
  cat ("\n\nCoefficiants:\n")
  print (object$coefs)
  cat ("\n\n")
}
  
