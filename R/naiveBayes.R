naiveBayes <- function(formula, data, ..., subset, na.action = na.pass) {
  call <- match.call()
  Yname <- as.character(formula[[2]])

  if (is.data.frame(data)) {
    ## handle formula
    m <- match.call(expand = FALSE)
    m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1)) 
      stop("naiveBayes cannot handle interaction terms")
    Y <- model.extract(m, "response")
    X <- m[,-attr(Terms, "response")]

    ## create tables
    apriori <- table(Y)
    tables <- lapply(X, function(x) {tab <- table(Y, x); tab / rowSums(tab)})

    ## fix dimname names
    for (i in 1:length(tables))
      names(dimnames(tables[[i]])) <- c(Yname, colnames(X)[i])
    names(dimnames(apriori)) <- Yname

    classlevels <- levels(Y)
  } else if (is.array(data)) {
    ## Find Class dimension
    Yind <- which(names(dimnames(data)) == Yname)

    ## Create Variable index
    deps <- strsplit(as.character(formula)[3], ".[+].")[[1]]
    if (length(deps) == 1 && deps == ".")
      deps <- names(dimnames(data))[-Yind]
    Vind <- which(names(dimnames(data)) %in% deps)
    
    ## create tables
    apriori <- margin.table(data, Yind)
    tables <- lapply(Vind,
                     function(i) margin.table(data, c(Yind, i)) / as.numeric(apriori))

    classlevels <- names(apriori)
  } else stop("naiveBayes handles data frames or arrays only")

  structure(list(apriori = apriori,
                 tables = tables,
                 levels = classlevels,
                 call   = call
                 ),
            
            class = "naiveBayes"
            )
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
                               ...) {
  type <- match.arg(type)
  nattribs <- ncol(newdata)
  L <- sapply(1:nrow(newdata), function(i) {
    ndata <- as.numeric(newdata[i,])
    L <- object$apriori *
      apply(sapply(1:nattribs, function(v) {
        nd <- ndata[v]
        if(is.na(nd))
          rep(1, length(object$apriori))
        else {
          prob <- object$tables[[v]][,nd]
          prob[prob == 0] <- threshold
          prob
        }
      }), 1, prod)
    L / sum(L)
  })
  if (type == "class")
    factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
  else
    t(L)
}
