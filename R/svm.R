svm <-
function (x,
          y,
          svm.type    = NULL,
          kernel.type = "radial",
          degree      = 3,
          gamma       = 1/ncol(as.matrix(x)),
          coef0       = 0,
          cost        = 1,
          nu          = 0.5,
          cachesize   = 40,
          tolerance   = 0.001,
          epsilon     = 0.5)
{
  if (is.vector(x))
    x <- t(t(x))
  else
    x <- as.matrix(x)

  if (is.null (svm.type)) svm.type <-
      if (is.factor(y)) "C-classification" else "regression"

  svm.type    <- pmatch (svm.type, c("C-classification",
                                     "nu-classification",
                                     "one-classification",
                                     "regression"),1) - 1
  
  kernel.type <- pmatch (kernel.type, c("linear",
                                        "polynomial",
                                        "radial",
                                        "sigmoid"),3) - 1

  if (!is.vector(y) && !is.factor (y)) stop ("y must be a vector or a factor.")
  if (length(y) != nrow(x)) stop ("x and y don't match.")

  if (cachesize < 0.1) cachesize <- 0.1
  
  lev <- NULL
  # in case of classification: map levels into {-1,1}
  if (is.factor(y)) {
    lev <- levels (y)
    y <- codes (y) * 2 - 3
  } else if (svm.type < 3) {
    lev <- levels (as.ordered (y))
    y <- codes (as.ordered(y)) * 2 - 3
  }

  if (length (lev) > 2) stop ("sorry, can't handle more than 2 classes !")
    
  cret <- .C ("svmtrain",
              as.double  (t(x)), as.integer(nrow (x)), as.integer(ncol(x)),
              as.double  (y),
              as.integer (svm.type),
              as.integer (kernel.type),
              as.double  (degree),
              as.double  (gamma),
              as.double  (coef0),
              as.double  (cost),
              as.double  (nu),
              as.double  (cachesize),
              as.double  (tolerance),
              as.double  (epsilon),
              nr    = integer (1),
              index = integer (nrow(x)),
              coefs = double  (nrow(x)),
              rho   = double  (1)
             )
  
  ret <- list (
               call        = match.call(),
               svm.type    = svm.type,
               kernel.type = kernel.type,
               cost        = cost,
               degree      = degree,
               gamma       = gamma,
               coef0       = coef0,
               nu          = nu,
               epsilon     = epsilon,
               
               levels      = lev,
               nr          = cret$nr,                  #number of sv
               sv          = t(t(x[cret$index==1,])),        #copy of sv
               index       = which (cret$index==1),    #indexes of sv in x
               rho         = cret$rho,                 #constant in decision function
               coefs       = cret$coefs[cret$index==1] #coefficiants of sv
              )
  class (ret) <- "svm"
  ret
} 

predict.svm <- function (model, x, type = "class") {

  type <- pmatch (type, c("class","raw"), 1)

  if (is.vector (x))
    x <- t(t(x))
  else
    x <- as.matrix(x)
    
  if (ncol(x) != ncol(model$sv)) stop ("test vector does not match model !")

  ret <- .C ("svmclassify",
             #model
             as.double  (t(model$sv)),
             as.integer (nrow(model$sv)), as.integer(ncol(model$sv)),
             as.double  (model$coefs),
             as.double  (model$rho),
             
             #parameter
             as.integer (model$svm.type),
             as.integer (model$kernel.type),
             as.double  (model$degree),
             as.double  (model$gamma),
             as.double  (model$coef0),

             #test matrix
             as.double (t(x)), as.integer (nrow(x)),
             
             #decision-values
             ret = double  (nrow(x))
            )$ret

  if ((type == 2) || is.null(model$levels))
    ret
  else {
    ret2 <- rep (model$levels[2],nrow(x))
    ret2 [which(ret < 0)] <- model$levels[1]
    ret2
  }
}

print.svm <- function (model) {
  cat ("\nCall:\n",deparse (model$call),"\n\n")
  cat ("Parameters:\n")
  cat ("   SVM-Type: ",c("C-classification",
                         "nu-classification",
                         "one-classification",
                         "regression")[model$svm.type+1],"\n")
  cat (" SVM-Kernel: ",c("linear",
                         "polynomial",
                         "radial",
                         "sigmoid")[model$kernel.type+1],"\n")
  cat ("       cost: ",model$cost,"\n")
  cat ("     degree: ",model$degree,"\n")
  cat ("      gamma: ",model$gamma,"\n")
  cat ("     coef.0: ",model$coef0,"\n")
  cat ("         nu: ",model$nu,"\n")
  cat ("    epsilon: ",model$epsilon,"\n")
  cat ("       cost: ",model$cost,"\n\n")
  
  cat ("\nNumber of Support Vectors:\n",model$nr,"\n\n")
  cat ("Rho:\n",model$rho,"\n\n")
  
}

summary.svm <- function (x) {
  print (x)
  cat ("Support Vectors:\n")
  print (x$sv)
  cat ("\n\nCoefficiants:\n")
  print (x$coefs)
  cat ("\n\n")
}
  
