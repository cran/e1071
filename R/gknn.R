gknn <- function(x, ...)
  UseMethod("gknn")

gknn.formula <-
function (formula, data = NULL, ..., subset, na.action = na.pass, scale = TRUE)
{
  call <- match.call()
  if (!inherits(formula, "formula"))
    stop("method is only for formula objects")
  m <- match.call(expand.dots = FALSE)
  if (inherits(eval.parent(m$data), "matrix"))
    m$data <- as.data.frame(eval.parent(m$data))
  m$... <- NULL
  m$scale <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m$na.action <- na.action
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, "response")
  attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")
  attr(x, "xlevels") <- .getXlevels(Terms, m)
  if (length(scale) == 1)
    scale <- rep(scale, ncol(x))
  if (any(scale)) {
    remove <- unique(c(which(labels(Terms) %in% names(attr(x, "contrasts"))),
                       which(!scale))
    )
    scale <- !attr(x, "assign") %in% remove
  }
  ret <- gknn.default (x, y, scale = scale, ..., na.action = na.action)
  ret$call <- call
  ret$call[[1]] <- as.name("gknn")
  ret$terms <- Terms
  ret$na.action <- attr(x, "na.action")
  class(ret) <- c("gknn.formula", class(ret))
  return (ret)
}


gknn.default <- function(x, y, k = 1, method = NULL,
                         scale = TRUE, use_all = TRUE,
                         FUN = mean,
                         ...)
{
    if (length(scale) == 1)
        scale <- rep(scale, ncol(x))
    if (is.numeric(x) && any(scale)) {
        tmp <- scale(x[,scale])
        x[,scale] <- tmp
        attr(x, "scaled:center") <- attr(tmp, "scaled:center")
        attr(x, "scaled:scale") <- attr(tmp, "scaled:scale")
    }
    structure(list(
              x = x,
              y = y,
              k = k,
              FUN = FUN,
              method = method,
              use_all = use_all,
              scaled = is.numeric(x) && any(scale),
              scale = scale),
        class = "gknn"
    )
}

predict.gknn <- function(object, newdata,
                         type = c("class", "votes", "prob"),
                         ...,
                         na.action = na.pass)
{
    if (missing(newdata))
        return(fitted(object))

    type = match.arg(type)

    if (inherits(object, "gknn.formula")) {
        if(is.null(colnames(newdata)))
            colnames(newdata) <- colnames(object$x)
        newdata <- na.action(newdata)
        act <- attr(newdata, "na.action")
        newdata <- model.matrix(delete.response(terms(object)),
                                as.data.frame(newdata),
                                xlev = attr(object$x, "xlevels"))
    } else {
        newdata <- na.action(as.matrix(newdata))
        act <- attr(newdata, "na.action")
    }

    if (object$scaled)
        newdata[,object$scale] <- scale(newdata[,object$scale, drop = FALSE],
                         center = attr(object$x, "scaled:center"),
                         scale = attr(object$x, "scaled:scale")
                         )
    d <- dist(object$x, newdata, method = object$method)
    FUN <- function(x) {
        o <- order(x)
        ks <- which(x[o][object$k] == x) ## check for ties on kth place
        if (!object$use_all) ks <- sample(c(ks, ks), 1) ## handle ties
        lab <- object$y[c(head(o[1:object$k], -1), ks)]
        if (is.numeric(lab))
            object$FUN(lab)
        else {
            tab <- table(lab)
            switch(type,
                   class = levels(object$y)[sample(rep(which(tab == max(tab)), 2), 1)], ## break class tie by random
                   prob = prop.table(tab),
                   tab)
        }
    }
    ret <- apply(d, 2, FUN)
    if (is.matrix(ret))
      t(ret)
    else
      if (is.numeric(object$y))
          napredict(act, ret)
       else
          napredict(act, factor(ret, levels = levels(object$y)))
}

print.gknn <- function(x, ...)
{
    cat("Object of class 'gknn'.\n")
}

fitted.gknn <- function(object, ...)
    napredict(object$na.action, object$y)
