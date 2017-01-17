naiveBayes <- function(x, ...)
    UseMethod("naiveBayes")

naiveBayes.default <- function(x, y, laplace = 0, ...) {
    call <- match.call()
    Yname <- deparse(substitute(y))
    x <- as.data.frame(x)

    ## estimation-function
    est <- function(var)
        if (is.numeric(var)) {
            cbind(tapply(var, y, mean, na.rm = TRUE),
                  tapply(var, y, sd, na.rm = TRUE))
        } else {
            tab <- table(y, var)
            (tab + laplace) / (rowSums(tab) + laplace * nlevels(var))
        }

    ## create tables
    apriori <- table(y)
    tables <- lapply(x, est)

    ## fix dimname names
    for (i in 1:length(tables))
        names(dimnames(tables[[i]])) <- c(Yname, colnames(x)[i])
    names(dimnames(apriori)) <- Yname

    structure(list(apriori = apriori,
                   tables = tables,
                   levels = levels(y),
                   call   = call
                   ),

              class = "naiveBayes"
              )
}

naiveBayes.formula <- function(formula, data, laplace = 0, ...,
                               subset, na.action = na.pass) {
    call <- match.call()
    Yname <- as.character(formula[[2]])

    if (is.data.frame(data)) {
        ## handle formula
        m <- match.call(expand.dots = FALSE)
        m$... <- NULL
        m$laplace = NULL
        m$na.action <- na.action
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        Terms <- attr(m, "terms")
        if (any(attr(Terms, "order") > 1))
            stop("naiveBayes cannot handle interaction terms")
        Y <- model.extract(m, "response")
        X <- m[,-attr(Terms, "response"), drop = FALSE]

        return(naiveBayes(X, Y, laplace = laplace, ...))
    } else if (is.array(data)) {
        nam <- names(dimnames(data))
        ## Find Class dimension
        Yind <- which(nam == Yname)

        ## Create Variable index
        deps <- strsplit(as.character(formula)[3], ".[+].")[[1]]
        if (length(deps) == 1 && deps == ".")
            deps <- nam[-Yind]
        Vind <- which(nam %in% deps)

        ## create tables
        apriori <- margin.table(data, Yind)
        tables <- lapply(Vind,
                         function(i) (margin.table(data, c(Yind, i)) + laplace) /
                         (as.numeric(apriori) + laplace * dim(data)[i]))
        names(tables) <- nam[Vind]

        structure(list(apriori = apriori,
                       tables = tables,
                       levels = names(apriori),
                       call   = call
                       ),

                  class = "naiveBayes"
                  )
    } else stop("naiveBayes formula interface handles data frames or arrays only")

}


print.naiveBayes <- function(x, ...) {
    cat("\nNaive Bayes Classifier for Discrete Predictors\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\nA-priori probabilities:\n")
    print(x$apriori / sum(x$apriori))

    cat("\nConditional probabilities:\n")
    for (i in x$tables) {print(i); cat("\n")}

}

predict.naiveBayes <- function(object,
                               newdata,
                               type = c("class", "raw"),
                               threshold = 0.001,
                               eps = 0,
                               ...) {
  # Initializing
  type <- match.arg(type)
  newdata <- as.data.frame(newdata)
  attribs <- match(names(object$tables), names(newdata))
  log_eps <- log(eps)
  logp <- lapply(object$tables, log)
  aprior <- log(object$apriori)

  # Function for assigning log probabilities to new data
  assign_probabilities <- function(v) {
    nb.id <- match(v, names(logp))
    p <- logp[[nb.id]][, newdata[, v]]
    p <- ifelse(p < log_eps, log(threshold), p)
    return(p)
  }

  # Actually assign log probabilities to new data
  a <- lapply(names(newdata)[attribs], assign_probabilities)

  # Compute likelihood by summing log probabilities (including prior)
  prior_p <- matrix(rep(aprior, ncol(a[[1]])), nrow(a[[1]]))
  posterior_p <- Reduce("+", a, prior_p)
  normalize_log_p <- function(x) sapply(x, function(y) {1 / sum(exp(x - y))})
  L <- apply(posterior_p, 2, normalize_log_p)

  # Return probabilities
  if (type == "class")
    factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
  else
    t(L)
}

