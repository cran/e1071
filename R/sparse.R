read.matrix.csr <- function (file, fac = TRUE, ncol = NULL) {
  library(methods)
  if(!require(SparseM))
    stop("Need `SparseM' package!")
  
  con <- file(file)
  open(con)
  y <- vector()
  x <- new("matrix.csr")
  x@ia <- 1
  i <- 1
  maxcol <- 1
  while (isOpen(con) & length(buf <- readLines (con, 1)) > 0) {
    s <- strsplit(buf, "[ ]+",extended=TRUE)[[1]]
    
    ## y
    if (length(grep(":", s[1])) == 0) {
      y[i] <- if (fac) s[1] else as.numeric(s[1])
      s <- s[-1]
    }

    ## x-values
    if (length(s)) {
      tmp <- do.call("rbind", strsplit(s, ":"))
      x@ra <- c(x@ra, as.numeric(tmp[,2]))
      x@ja <- c(x@ja, as.numeric(tmp[,1]))
    }
    i <- i + 1
    x@ia[i] <- x@ia[i-1] + length(s) 
  }
  x@dimension <- c(i - 1, if (is.null(ncol)) max(x@ja) else max(ncol, x@ja))
  class(x) <- "matrix.csr"
  
  if (length(y)) 
    list (x = x, y = if (fac) as.factor(y) else y)
  else x
}

write.matrix.csr <- function (x, file="out.dat", y=NULL) {
  library(methods)
  if (!is.null(y) & (length(y) != nrow(x)))
    stop(paste("Length of y (=", length(y),
                 ") does not match number of rows of x (=",
                 nrow(x), ")!", sep=""))
  sink(file)
  for (i in 1:nrow(x)) {
    if (!is.null(y)) cat (y[i],"")
    for (j in x@ia[i]:(x@ia[i+1] - 1))
      cat(x@ja[j], ":", x@ra[j], " ", sep="")
    cat("\n")
  }
  sink()
}

na.fail.matrix.csr <- function(object, ...) {
  library(methods)
  if (any(is.na(object@ra)))
    stop("missing values in object") else return(object)
}






